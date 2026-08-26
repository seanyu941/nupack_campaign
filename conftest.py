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
TARGET_SEQUENCE = "AAGCGCGGAAGC"

# (trunc_5, trunc_3, mutations) where mutations is a list of (pos, from, to, type)
FIXTURE_ROWS = [
    (0, 0, []),
    (0, 0, [(2, "T", "A", "transversion"), (5, "C", "G", "transversion")]),
    (0, 0, [(3, "T", "C", "transition"), (7, "G", "A", "transition")]),
    (1, 0, []),
    (1, 0, [(2, "T", "A", "transversion"), (5, "C", "G", "transversion")]),
    (1, 0, [(3, "T", "C", "transition"), (7, "G", "A", "transition")]),
    (1, 1, []),
    (1, 1, [(2, "T", "A", "transversion"), (5, "C", "G", "transversion")]),
    (1, 1, [(3, "T", "C", "transition"), (7, "G", "A", "transition")]),
]

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

        for slot, (position, original, mutated, kind) in enumerate(mutations, start=1):
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
            num_transitions=sum(m[3] == "transition" for m in mutations),
            num_transversions=sum(m[3] == "transversion" for m in mutations),
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
def thermo_engine() -> StubEngine:
    return StubEngine(seed="test")


@pytest.fixture
def grid() -> dict:
    return {key: list(values) for key, values in FIXTURE_GRID.items()}
