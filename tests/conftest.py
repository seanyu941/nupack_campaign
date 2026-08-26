"""Shared fixtures.

Everything runs against an in-memory SQLite database with a small scan, so the
suite finishes in about a second and leaves nothing behind. The schema is the
same one the real database uses, since both come from the ORM metadata.

The fixture scan is a miniature of the real thing: one 12 nt parent, three
truncation combinations, and an unmutated reference plus two double mutants at
each. That is enough for the per-truncation ddG join and the per-class cap to
have something to work on.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from sqlalchemy.orm import Session

from nupack_campaign.catalog import classify, gc_content, mutation_signature
from nupack_campaign.db import init_db, make_engine, make_session_factory
from nupack_campaign.engines import StubEngine
from nupack_campaign.models import Target, Variant

REPO_ROOT = Path(__file__).resolve().parents[1]

PARENT = "GCTTCCGCGCTT"
# The parent's reverse complement (AAGCGCGGAAGC) sits at offset 6, with flanks
# either side, so the aligner has to find the window rather than assume the two
# sequences start together. That mirrors the real 76 nt target, where the 25 nt
# design binds at offset 51.
TARGET_SEQUENCE = "ACAATC" + "AAGCGCGGAAGC" + "TCATAC"
BINDING_OFFSET = 6

# (trunc_5prime, trunc_3prime, [(position, substitution kind), ...]).
#
# Only the position and the kind are fixed here. The original base is read off
# the truncated sequence and the substituted base derived from it, so a mutation
# always actually changes the sequence. Hardcoding the bases instead lets a
# "mutation" silently set a base to the value it already had once the sequence
# is truncated, which makes the paired-base counts wrong in a way that is hard
# to spot.
FIXTURE_ROWS = [
    (0, 0, []),
    (0, 0, [(2, "transversion"), (5, "transversion")]),
    (0, 0, [(3, "transition"), (7, "transition")]),
    (1, 0, []),
    (1, 0, [(2, "transversion"), (5, "transversion")]),
    (1, 0, [(3, "transition"), (7, "transition")]),
    (1, 1, []),
    (1, 1, [(2, "transversion"), (5, "transversion")]),
    (1, 1, [(3, "transition"), (7, "transition")]),
]

TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
TRANSVERSION = {"A": "C", "G": "T", "C": "A", "T": "G"}


def substitute(base: str, kind: str) -> str:
    return TRANSITION[base] if kind == "transition" else TRANSVERSION[base]


FIXTURE_GRID = {
    "temperature_c": [25.0, 37.0],
    "na_molar": [0.15],
    "mg_molar": [0.0],
    "material": ["dna04"],
    "ensemble": ["stacking"],
}


def build_fixture_variants() -> list[Variant]:
    variants = []
    for trunc_5, trunc_3, mutations in FIXTURE_ROWS:
        sequence = list(PARENT[trunc_5 : len(PARENT) - trunc_3 if trunc_3 else None])
        record: dict = {}

        for slot, (position, kind) in enumerate(mutations, start=1):
            original = sequence[position]
            mutated = substitute(original, kind)
            sequence[position] = mutated
            record[f"position_{slot}"] = position
            record[f"original_base_{slot}"] = original
            record[f"mutated_base_{slot}"] = mutated
            record[f"mutation_type_{slot}"] = kind

        for slot in (1, 2):
            record.setdefault(f"position_{slot}", None)
            record.setdefault(f"original_base_{slot}", None)
            record.setdefault(f"mutated_base_{slot}", None)
            record.setdefault(f"mutation_type_{slot}", "none")

        joined = "".join(sequence)
        trunc_total = trunc_5 + trunc_3
        record.update(
            sequence=joined,
            length_nt=len(joined),
            gc_content=gc_content(joined),
            trunc_5prime=trunc_5,
            trunc_3prime=trunc_3,
            trunc_total=trunc_total,
            num_mutations=len(mutations),
            num_transitions=sum(kind == "transition" for _, kind in mutations),
            num_transversions=sum(kind == "transversion" for _, kind in mutations),
            variant_class=classify(len(mutations), trunc_total),
            is_reference=int(len(mutations) == 0 and trunc_total == 0),
        )
        record["mutation_signature"] = mutation_signature(record)
        record["name"] = f"t{trunc_5}_{trunc_3}_" + (record["mutation_signature"] or "wt")
        variants.append(Variant(**record))

    return variants


@pytest.fixture
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture
def query_path(repo_root: Path) -> Path:
    return repo_root / "sql" / "select_candidates.sql"


@pytest.fixture
def strongest_query_path(repo_root: Path) -> Path:
    return repo_root / "sql" / "select_strongest_binders.sql"


@pytest.fixture
def trend_query_path(repo_root: Path) -> Path:
    return repo_root / "sql" / "trends" / "trend_by_dimension.sql"


@pytest.fixture
def engine():
    db_engine = make_engine("sqlite:///:memory:")
    init_db(db_engine)
    return db_engine


@pytest.fixture
def session(engine) -> Session:
    factory = make_session_factory(engine)
    with factory() as db_session:
        yield db_session


@pytest.fixture
def populated_session(session: Session) -> Session:
    session.add_all(build_fixture_variants())
    session.commit()
    return session


@pytest.fixture
def target(populated_session: Session) -> Target:
    record = Target(name="fixture_target", sequence=TARGET_SEQUENCE)
    populated_session.add(record)
    populated_session.commit()
    return record


@pytest.fixture
def aligned_target(populated_session: Session, target: Target) -> Target:
    from nupack_campaign.catalog import align_variants

    align_variants(populated_session, target)
    return target


@pytest.fixture
def thermo_engine() -> StubEngine:
    return StubEngine(seed="test")


@pytest.fixture
def grid() -> dict:
    return {key: list(values) for key, values in FIXTURE_GRID.items()}


def fixture_name(trunc_5prime: int, trunc_3prime: int, mutation_index: int = 0) -> str:
    """Look up a fixture variant by its scan coordinates.

    Names carry the substituted bases, which are derived rather than fixed, so
    tests ask for "the second mutant at this truncation" instead of hardcoding a
    signature that would change if the parent sequence changed.
    """
    matching = [
        variant.name
        for variant in build_fixture_variants()
        if variant.trunc_5prime == trunc_5prime and variant.trunc_3prime == trunc_3prime
    ]
    return matching[mutation_index]
