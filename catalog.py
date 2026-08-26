"""Load a design scan CSV into the ``variants`` table.

The expected columns are the ones a mutation and truncation scan produces:

    sequence, length, gc_content, trunc_5prime, trunc_3prime,
    position_1, original_base_1, mutated_base_1, mutation_type_1,
    position_2, original_base_2, mutated_base_2, mutation_type_2,
    num_transitions, num_transversions

Only ``sequence`` is required. Everything else is filled in from the sequence or
defaulted, so a plain one-column list of sequences also loads. Any free energy
columns in the file are ignored here, since those are results rather than design
metadata; ``importer.py`` is what reads those.

Four columns are derived rather than read, so they cannot disagree with the row
they describe: ``length_nt``, ``trunc_total``, ``num_mutations`` and
``variant_class``. ``gc_content`` is recomputed from the sequence and checked
against the file if the file supplies it.

Loading is idempotent. Variants already present by name are left alone, and a
name that comes back with a different sequence is rejected, since that would
silently invalidate every result already tied to that variant_id.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import Session

from .models import Variant

VALID_BASES = set("ACGT")
GC_TOLERANCE = 0.01

SCAN_COLUMNS = (
    "length",
    "gc_content",
    "trunc_5prime",
    "trunc_3prime",
    "position_1",
    "original_base_1",
    "mutated_base_1",
    "mutation_type_1",
    "position_2",
    "original_base_2",
    "mutated_base_2",
    "mutation_type_2",
    "num_transitions",
    "num_transversions",
)


class CatalogError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class LoadStats:
    inserted: int
    skipped: int
    by_class: dict[str, int]

    def __str__(self) -> str:
        classes = ", ".join(f"{name} {count}" for name, count in sorted(self.by_class.items()))
        return f"{self.inserted} inserted, {self.skipped} already present ({classes})"


def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(base in "GC" for base in sequence.upper()) / len(sequence)


def _clean_optional(value) -> str | None:
    """Blank, NaN and the literal string 'none' all mean "no mutation here"."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "na"}:
        return None
    return text


def _clean_position(value) -> int | None:
    text = _clean_optional(value)
    if text is None:
        return None
    return int(float(text))


def classify(num_mutations: int, trunc_total: int) -> str:
    """Bucket a row for grouping and for the per-class cap in selection.

    Mutation count and truncation are independent axes of the scan, so the class
    names carry both. Trend queries that care about only one axis group by
    ``num_mutations`` or ``trunc_total`` directly instead.
    """
    if num_mutations == 0:
        return "reference" if trunc_total == 0 else "truncation_only"
    stem = {1: "single_mutation", 2: "double_mutation"}[num_mutations]
    return stem if trunc_total == 0 else f"{stem}_truncated"


def mutation_signature(row: dict) -> str:
    """Readable label such as "G0A,C1A", empty for an unmutated row."""
    parts = []
    for slot in (1, 2):
        position = row.get(f"position_{slot}")
        original = row.get(f"original_base_{slot}")
        mutated = row.get(f"mutated_base_{slot}")
        if position is not None and original and mutated:
            parts.append(f"{original}{position}{mutated}")
    return ",".join(parts)


def read_catalog(path: str | Path) -> pd.DataFrame:
    """Read the CSV and normalise it into the columns the loader writes."""
    frame = pd.read_csv(path)

    if "sequence" not in frame.columns:
        raise CatalogError(f"{path} has no 'sequence' column")

    for column in SCAN_COLUMNS:
        if column not in frame.columns:
            frame[column] = None

    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()

    bad = frame.loc[~frame["sequence"].apply(lambda s: bool(s) and set(s) <= VALID_BASES)]
    if len(bad):
        raise CatalogError(
            f"non-ACGT characters in sequences at rows {bad.index[:5].tolist()}. "
            "Degenerate bases are not supported."
        )

    records = []
    for index, raw in enumerate(frame.to_dict("records")):
        sequence = raw["sequence"]

        record = {
            "sequence": sequence,
            "length_nt": len(sequence),
            "gc_content": gc_content(sequence),
            "trunc_5prime": int(raw["trunc_5prime"] or 0),
            "trunc_3prime": int(raw["trunc_3prime"] or 0),
            "num_transitions": int(raw["num_transitions"] or 0),
            "num_transversions": int(raw["num_transversions"] or 0),
        }

        for slot in (1, 2):
            record[f"position_{slot}"] = _clean_position(raw[f"position_{slot}"])
            record[f"original_base_{slot}"] = _clean_optional(raw[f"original_base_{slot}"])
            record[f"mutated_base_{slot}"] = _clean_optional(raw[f"mutated_base_{slot}"])
            record[f"mutation_type_{slot}"] = (
                _clean_optional(raw[f"mutation_type_{slot}"]) or "none"
            )

        # A slot counts as mutated when it names a position and a substitution,
        # not when mutation_type happens to be filled in.
        record["num_mutations"] = sum(
            record[f"position_{slot}"] is not None and record[f"mutated_base_{slot}"] is not None
            for slot in (1, 2)
        )
        record["trunc_total"] = record["trunc_5prime"] + record["trunc_3prime"]
        record["variant_class"] = classify(record["num_mutations"], record["trunc_total"])
        record["mutation_signature"] = mutation_signature(record)
        record["is_reference"] = int(
            record["num_mutations"] == 0 and record["trunc_total"] == 0
        )

        # The file's own length and gc_content are treated as a cross check
        # rather than as input, so a stale column cannot poison a trend query.
        if raw["length"] is not None and not pd.isna(raw["length"]):
            if int(raw["length"]) != record["length_nt"]:
                raise CatalogError(
                    f"row {index}: length column says {int(raw['length'])} but the "
                    f"sequence is {record['length_nt']} nt"
                )
        if raw["gc_content"] is not None and not pd.isna(raw["gc_content"]):
            if abs(float(raw["gc_content"]) - record["gc_content"]) > GC_TOLERANCE:
                raise CatalogError(
                    f"row {index}: gc_content column says {float(raw['gc_content']):.4f} "
                    f"but the sequence gives {record['gc_content']:.4f}"
                )

        records.append(record)

    catalog = pd.DataFrame(records)

    if "name" in frame.columns:
        catalog["name"] = frame["name"].astype(str).str.strip().to_numpy()
        duplicates = catalog["name"][catalog["name"].duplicated()]
        if len(duplicates):
            raise CatalogError(
                f"duplicate names in {path}: {duplicates.head(3).tolist()}"
            )
    else:
        catalog["name"] = build_names(catalog)

    duplicates = catalog["name"][catalog["name"].duplicated()]
    if len(duplicates):
        raise CatalogError(
            f"the scan metadata does not uniquely identify these rows: "
            f"{duplicates.head(3).tolist()}. Add a 'name' column to the CSV."
        )

    return catalog


def build_names(catalog: pd.DataFrame) -> pd.Series:
    """Build a stable name from the scan coordinates.

    Something like ``t0_0_wt`` or ``t2_1_G0A,C1A``. Deriving the name from the
    scan metadata rather than the row order means reloading a re-sorted file
    matches the existing variants instead of duplicating them.

    A file with no scan metadata at all, such as a plain list of sequences,
    would give every row the same name. When that happens a short sequence hash
    is appended to every row rather than only to the colliding ones, so the
    names stay stable if the file later gains another row.
    """
    trunc = (
        "t"
        + catalog["trunc_5prime"].astype(str)
        + "_"
        + catalog["trunc_3prime"].astype(str)
    )
    names = trunc + "_" + catalog["mutation_signature"].replace("", "wt")

    if names.duplicated().any():
        digest = catalog["sequence"].map(
            lambda s: hashlib.blake2b(s.encode(), digest_size=3).hexdigest()
        )
        names = names + "_" + digest

    return names


def load_variants(session: Session, path: str | Path, batch_size: int = 5000) -> LoadStats:
    catalog = read_catalog(path)

    existing = dict(session.execute(select(Variant.name, Variant.sequence)).all())

    pending: list[Variant] = []
    inserted = 0
    skipped = 0
    by_class: dict[str, int] = {}

    for row in catalog.to_dict("records"):
        name = row["name"]

        if name in existing:
            if existing[name] != row["sequence"]:
                raise CatalogError(
                    f"variant {name} is already in the database with a different "
                    "sequence. Give the new design its own name so existing results "
                    "keep pointing at the sequence they were computed for."
                )
            skipped += 1
            continue

        pending.append(Variant(**row))
        inserted += 1
        by_class[row["variant_class"]] = by_class.get(row["variant_class"], 0) + 1

        if len(pending) >= batch_size:
            session.add_all(pending)
            session.flush()
            pending.clear()

    if pending:
        session.add_all(pending)

    session.commit()
    return LoadStats(inserted=inserted, skipped=skipped, by_class=by_class)
