"""Import a scan CSV that already has free energies in it.

A scan CSV stores results wide, one group of columns per temperature:

    cdna_dg_25C, complex_dg_25C, delta_g_binding_25C,
    cdna_dg_37C, complex_dg_37C, delta_g_binding_37C, ...

That shape is fine to look at and awkward to query. Adding a temperature means
adding three more columns and rewriting every query that mentions them, and
asking "how does binding change with temperature" means unpivoting by hand
every time. This module unpivots the file into ``results``, one row per
(variant, condition), which is the shape the rest of the pipeline expects.

It also recovers the target strand energy, which the file does not carry. Since

    delta_g_binding = complex_dg - cdna_dg - target_dg

the target term is whatever the other three imply, and it should come out the
same for every row at a given temperature because it does not depend on the
variant. The importer checks that it does, within a tolerance, and refuses the
file if it does not. That doubles as a check that the three columns are
self-consistent.

Importing creates a run like any other, so imported and simulated results sit
side by side and every downstream query works on both.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import Session

from .catalog import align_variants
from .db import analyze_database
from .engines import check_binding_site
from .models import Condition, Result, Run, Target, TargetEnergy, Variant
from .sweep import current_git_sha

# Matches cdna_dg_25C, complex_dg_37C, delta_g_binding_25C and so on.
COLUMN_PATTERN = re.compile(
    r"^(?P<kind>cdna_dg|complex_dg|delta_g_binding)_(?P<temperature>-?\d+(?:\.\d+)?)C$"
)

# The scan CSV rounds to 3 decimals, so the implied target energy wobbles in the
# 4th. Anything larger means the three columns disagree with each other.
TARGET_ENERGY_TOLERANCE = 0.01


class ImportError_(ValueError):
    """Raised when the file cannot be reconciled with the schema."""


@dataclass(frozen=True, slots=True)
class ImportStats:
    run_id: int
    n_rows: int
    n_conditions: int
    temperatures: tuple[float, ...]
    target_energies: dict[float, float]
    binding_offset: int = 0
    binding_paired: int = 0

    def __str__(self) -> str:
        temps = ", ".join(f"{t:g}C" for t in self.temperatures)
        return (
            f"run {self.run_id}: {self.n_rows} rows across "
            f"{self.n_conditions} conditions ({temps})"
        )


def find_energy_columns(columns: list[str]) -> dict[float, dict[str, str]]:
    """Group the wide energy columns by the temperature in their name."""
    found: dict[float, dict[str, str]] = {}
    for column in columns:
        match = COLUMN_PATTERN.match(column)
        if match:
            temperature = float(match.group("temperature"))
            found.setdefault(temperature, {})[match.group("kind")] = column

    complete = {}
    for temperature, kinds in sorted(found.items()):
        missing = {"cdna_dg", "complex_dg", "delta_g_binding"} - set(kinds)
        if missing:
            raise ImportError_(
                f"{temperature:g}C is missing {', '.join(sorted(missing))}. "
                "All three energy columns are needed to recover the target term."
            )
        complete[temperature] = kinds

    if not complete:
        raise ImportError_(
            "no energy columns found. Expected names like cdna_dg_37C, "
            "complex_dg_37C and delta_g_binding_37C."
        )
    return complete


def recover_target_energy(frame: pd.DataFrame, columns: dict[str, str]) -> float:
    """Back out the target strand energy and check it is actually constant."""
    implied = (
        frame[columns["complex_dg"]]
        - frame[columns["cdna_dg"]]
        - frame[columns["delta_g_binding"]]
    )
    spread = float(implied.max() - implied.min())

    if spread > TARGET_ENERGY_TOLERANCE:
        raise ImportError_(
            f"the implied target energy varies by {spread:.4f} kcal/mol across the "
            f"file, which it should not since it does not depend on the variant. "
            f"Either the three energy columns disagree or the file mixes targets."
        )

    return float(implied.mean())


def import_scan_csv(
    session: Session,
    path: str | Path,
    target_name: str,
    target_sequence: str,
    na_molar: float = 0.15,
    mg_molar: float = 0.0,
    material: str = "dna04",
    ensemble: str = "stacking",
    note: str | None = None,
    batch_size: int = 5000,
) -> ImportStats:
    """Unpivot a scan CSV into ``results`` under a new run.

    Variants must already be loaded. The importer matches rows to variants by
    sequence, so the same file can be used for both stages.
    """
    frame = pd.read_csv(path)
    energy_columns = find_energy_columns(list(frame.columns))

    if "sequence" not in frame.columns:
        raise ImportError_(f"{path} has no 'sequence' column to match variants on")
    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()

    variant_ids = dict(session.execute(select(Variant.sequence, Variant.variant_id)).all())
    unknown = set(frame["sequence"]) - set(variant_ids)
    if unknown:
        raise ImportError_(
            f"{len(unknown)} sequences in the file are not in the variants table. "
            "Run load-variants on this file first."
        )

    # Confirm the designs have somewhere to bind before importing anything.
    # A wrong --target-sequence would otherwise produce a full set of results
    # that look fine and mean nothing.
    site = check_binding_site(list(frame["sequence"].unique()), target_sequence)

    target = session.execute(
        select(Target).where(Target.name == target_name)
    ).scalar_one_or_none()
    if target is None:
        target = Target(name=target_name, sequence=target_sequence.strip().upper())
        session.add(target)
        session.flush()

    align_variants(session, target)

    run = Run(
        git_sha=current_git_sha(),
        engine_name="imported",
        engine_version=Path(path).name,
        grid_hash="imported",
        status="running",
        note=note or f"imported from {Path(path).name}",
    )
    session.add(run)
    session.flush()

    target_energies: dict[float, float] = {}
    total_rows = 0

    for temperature, columns in energy_columns.items():
        target_dg = recover_target_energy(frame, columns)
        target_energies[temperature] = target_dg

        condition = _get_or_create_condition(
            session, temperature, na_molar, mg_molar, material, ensemble
        )

        session.add(
            TargetEnergy(
                run_id=run.run_id,
                target_id=target.target_id,
                condition_id=condition.condition_id,
                free_energy_kcal=target_dg,
            )
        )

        pending: list[Result] = []
        for row in frame[
            ["sequence", columns["cdna_dg"], columns["complex_dg"], columns["delta_g_binding"]]
        ].itertuples(index=False, name=None):
            sequence, cdna_dg, complex_dg, delta_g_binding = row
            pending.append(
                Result(
                    run_id=run.run_id,
                    variant_id=variant_ids[sequence],
                    target_id=target.target_id,
                    condition_id=condition.condition_id,
                    cdna_dg_kcal=float(cdna_dg),
                    complex_dg_kcal=float(complex_dg),
                    delta_g_binding_kcal=float(delta_g_binding),
                )
            )
            if len(pending) >= batch_size:
                session.add_all(pending)
                session.flush()
                total_rows += len(pending)
                pending = []

        if pending:
            session.add_all(pending)
            session.flush()
            total_rows += len(pending)

    run.status = "complete"
    from .models import utcnow

    run.finished_at = utcnow()
    session.commit()

    # Bulk load done, so refresh planner stats before anyone queries this.
    analyze_database(session)

    return ImportStats(
        run_id=run.run_id,
        n_rows=total_rows,
        n_conditions=len(energy_columns),
        temperatures=tuple(sorted(energy_columns)),
        target_energies=target_energies,
        binding_offset=site.offset,
        binding_paired=site.n_complementary,
    )


def _get_or_create_condition(
    session: Session,
    temperature_c: float,
    na_molar: float,
    mg_molar: float,
    material: str,
    ensemble: str,
) -> Condition:
    condition = session.execute(
        select(Condition).where(
            Condition.temperature_c == temperature_c,
            Condition.na_molar == na_molar,
            Condition.mg_molar == mg_molar,
            Condition.material == material,
            Condition.ensemble == ensemble,
        )
    ).scalar_one_or_none()

    if condition is None:
        condition = Condition(
            temperature_c=temperature_c,
            na_molar=na_molar,
            mg_molar=mg_molar,
            material=material,
            ensemble=ensemble,
        )
        session.add(condition)
        session.flush()

    return condition
