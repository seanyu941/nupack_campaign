"""Parameter sweep: expand the grid, work out what is missing, fill it in.

The grid lives in ``config/sweep_grid.yaml`` rather than in code, so widening a
sweep is a config change and the file can be hashed and stored against the run.

The runner is resumable. It never asks "have I already done this?" in Python.
It asks the database for the cells that have no row yet and computes those, so
an interrupted sweep picks up where it stopped and a finished sweep re-run is a
no-op.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
import time
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml
from sqlalchemy import and_, func, literal, select
from sqlalchemy.orm import Session

from .catalog import align_variants
from .db import analyze_database
from .engines import ConditionPoint, ThermoEngine, check_binding_site
from .models import Condition, Result, Run, Target, TargetEnergy, Variant, utcnow

GRID_KEYS = ("temperature_c", "na_molar", "mg_molar", "material", "ensemble")


@dataclass(frozen=True, slots=True)
class SweepStats:
    run_id: int
    total_cells: int
    computed_cells: int
    skipped_cells: int
    elapsed_seconds: float

    @property
    def cells_per_second(self) -> float:
        if self.elapsed_seconds <= 0:
            return 0.0
        return self.computed_cells / self.elapsed_seconds


def load_grid_file(path: str | Path) -> dict[str, list]:
    """Read the grid file and check it has the keys the model needs."""
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict):
        raise ValueError(f"{path} should be a mapping of parameter name to list of values")

    missing = [key for key in GRID_KEYS if key not in raw]
    if missing:
        raise ValueError(f"{path} is missing required keys: {', '.join(missing)}")

    grid: dict[str, list] = {}
    for key in GRID_KEYS:
        values = raw[key]
        if not isinstance(values, list) or not values:
            raise ValueError(f"{path}: {key} should be a non-empty list")
        # Dedupe while keeping the order written in the file.
        grid[key] = list(dict.fromkeys(values))
    return grid


def expand_grid(grid: dict[str, list]) -> list[ConditionPoint]:
    """Full cross product of the grid, in a stable order."""
    combinations = itertools.product(*(grid[key] for key in GRID_KEYS))
    return [
        ConditionPoint(
            temperature_c=float(temperature),
            na_molar=float(na),
            mg_molar=float(mg),
            material=str(material),
            ensemble=str(ensemble),
        )
        for temperature, na, mg, material, ensemble in combinations
    ]


def grid_hash(grid: dict[str, list]) -> str:
    """Stable hash of the grid contents, stored on the run row."""
    canonical = json.dumps(grid, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def current_git_sha(default: str = "unknown") -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return default


def ensure_conditions(
    session: Session, points: Sequence[ConditionPoint]
) -> dict[tuple, int]:
    """Insert any grid points that are not in ``conditions`` yet.

    Returns a mapping from condition key to condition_id covering every point
    passed in, whether it was already there or was just created.
    """
    existing = {
        (row.temperature_c, row.na_molar, row.mg_molar, row.material, row.ensemble):
            row.condition_id
        for row in session.execute(select(Condition)).scalars()
    }

    new_rows = []
    for point in points:
        if point.as_key() not in existing:
            new_rows.append(
                Condition(
                    temperature_c=point.temperature_c,
                    na_molar=point.na_molar,
                    mg_molar=point.mg_molar,
                    material=point.material,
                    ensemble=point.ensemble,
                )
            )

    if new_rows:
        session.add_all(new_rows)
        session.flush()
        for row in new_rows:
            key = (row.temperature_c, row.na_molar, row.mg_molar, row.material, row.ensemble)
            existing[key] = row.condition_id

    return {point.as_key(): existing[point.as_key()] for point in points}


def start_run(
    session: Session,
    engine: ThermoEngine,
    grid: dict[str, list],
    note: str | None = None,
) -> Run:
    run = Run(
        git_sha=current_git_sha(),
        engine_name=engine.name,
        engine_version=engine.version,
        grid_hash=grid_hash(grid),
        status="running",
        note=note,
    )
    session.add(run)
    session.flush()
    return run


def resume_run(session: Session, run_id: int, grid: dict[str, list]) -> Run:
    """Reuse an existing run row, refusing if the grid file has moved on."""
    run = session.get(Run, run_id)
    if run is None:
        raise ValueError(f"No run with id {run_id}")
    if run.grid_hash != grid_hash(grid):
        raise ValueError(
            f"Run {run_id} was started against a different sweep grid. Start a new "
            f"run instead so the two grids do not end up mixed in one result set."
        )
    run.status = "running"
    run.finished_at = None
    return run


def pending_pairs(
    session: Session,
    run_id: int,
    target_id: int,
    condition_ids: Sequence[int] | None = None,
) -> list[tuple[int, int]]:
    """Every (variant_id, condition_id) with no row in ``results`` for this run.

    This is a cross join of the two catalogs against an anti-join on results.
    Doing it in SQL keeps the whole "what is left to do" question on one side of
    the wire, and it stays cheap as the sweep grows because the outer join hits
    uq_results_cell rather than scanning.
    """
    stmt = (
        select(Variant.variant_id, Condition.condition_id)
        .select_from(Variant)
        .join(Condition, literal(True))
        .outerjoin(
            Result,
            and_(
                Result.run_id == run_id,
                Result.target_id == target_id,
                Result.variant_id == Variant.variant_id,
                Result.condition_id == Condition.condition_id,
            ),
        )
        .where(Result.result_id.is_(None))
        .order_by(Variant.variant_id, Condition.condition_id)
    )
    if condition_ids is not None:
        stmt = stmt.where(Condition.condition_id.in_(list(condition_ids)))

    return [tuple(row) for row in session.execute(stmt).all()]


def _chunked(items: Iterable, size: int) -> Iterator[list]:
    batch: list = []
    for item in items:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def ensure_target_energies(
    session: Session,
    thermo_engine: ThermoEngine,
    run: Run,
    target: Target,
    key_to_condition_id: dict[tuple, int],
    points: Sequence[ConditionPoint],
) -> dict[int, float]:
    """Evaluate the target strand once per condition.

    The target term in delta_g_binding does not depend on the variant, so this
    is one evaluation per condition rather than one per cell. On a 54,660 row
    scan across two conditions that is 2 evaluations instead of 109,320.
    """
    existing = {
        row.condition_id: row.free_energy_kcal
        for row in session.execute(
            select(TargetEnergy).where(
                TargetEnergy.run_id == run.run_id, TargetEnergy.target_id == target.target_id
            )
        ).scalars()
    }

    for point in points:
        condition_id = key_to_condition_id[point.as_key()]
        if condition_id in existing:
            continue
        free_energy = thermo_engine.evaluate_strand(target.sequence, point)
        session.add(
            TargetEnergy(
                run_id=run.run_id,
                target_id=target.target_id,
                condition_id=condition_id,
                free_energy_kcal=free_energy,
            )
        )
        existing[condition_id] = free_energy

    session.commit()
    return existing


def run_sweep(
    session: Session,
    thermo_engine: ThermoEngine,
    run: Run,
    target: Target,
    points: Sequence[ConditionPoint],
    batch_size: int = 500,
    progress_every: int = 5000,
    verbose: bool = True,
) -> SweepStats:
    """Fill in the missing cells for ``run``.

    Rows are committed in batches. Per-row commits dominate the runtime once the
    grid gets past a few thousand cells, and the stub engine is fast enough that
    the commit is the expensive part.
    """
    key_to_condition_id = ensure_conditions(session, points)
    condition_ids = list(key_to_condition_id.values())
    session.commit()

    variant_sequences_all = [
        row[0] for row in session.execute(select(Variant.sequence)).all()
    ]
    if variant_sequences_all:
        check_binding_site(variant_sequences_all, target.sequence)
        align_variants(session, target)

    target_dg = ensure_target_energies(
        session, thermo_engine, run, target, key_to_condition_id, points
    )

    variant_sequences = dict(
        session.execute(select(Variant.variant_id, Variant.sequence)).all()
    )
    if not variant_sequences:
        raise RuntimeError(
            "No variants in the database. Load the catalog first with "
            "`python -m nupack_campaign.cli load-variants`."
        )

    condition_by_id = {
        condition_id: point for point, condition_id in
        ((p, key_to_condition_id[p.as_key()]) for p in points)
    }

    total_cells = len(variant_sequences) * len(condition_ids)
    todo = pending_pairs(session, run.run_id, target.target_id, condition_ids)
    skipped = total_cells - len(todo)

    if verbose:
        print(
            f"[sweep] run {run.run_id}: {len(variant_sequences)} variants x "
            f"{len(condition_ids)} conditions = {total_cells} cells, "
            f"{len(todo)} to compute, {skipped} already present"
        )

    started = time.perf_counter()
    done = 0

    for batch in _chunked(todo, batch_size):
        rows = []
        for variant_id, condition_id in batch:
            point = condition_by_id[condition_id]
            cell_started = time.perf_counter()
            result = thermo_engine.evaluate_binding(
                variant_sequences[variant_id],
                target.sequence,
                point,
                target_dg[condition_id],
            )
            elapsed = time.perf_counter() - cell_started

            rows.append(
                Result(
                    run_id=run.run_id,
                    variant_id=variant_id,
                    target_id=target.target_id,
                    condition_id=condition_id,
                    cdna_dg_kcal=result.cdna_dg_kcal,
                    complex_dg_kcal=result.complex_dg_kcal,
                    delta_g_binding_kcal=result.delta_g_binding_kcal,
                    mfe_kcal=result.mfe_kcal,
                    ensemble_defect=result.ensemble_defect,
                    compute_seconds=elapsed,
                )
            )

        session.add_all(rows)
        session.commit()
        done += len(rows)

        if verbose and progress_every and done % progress_every < batch_size:
            pct = 100.0 * done / max(len(todo), 1)
            print(f"[sweep] {done}/{len(todo)} cells ({pct:.1f}%)")

    elapsed_total = time.perf_counter() - started
    run.status = "complete"
    run.finished_at = utcnow()
    session.commit()

    if done:
        analyze_database(session)

    if verbose:
        print(f"[sweep] run {run.run_id} complete in {elapsed_total:.2f}s")

    return SweepStats(
        run_id=run.run_id,
        total_cells=total_cells,
        computed_cells=done,
        skipped_cells=skipped,
        elapsed_seconds=elapsed_total,
    )


def run_summary(session: Session, run_id: int) -> dict:
    """Small aggregate used by the CLI and the notebook."""
    row = session.execute(
        select(
            func.count(Result.result_id),
            func.count(func.distinct(Result.variant_id)),
            func.count(func.distinct(Result.condition_id)),
            func.min(Result.delta_g_binding_kcal),
            func.max(Result.delta_g_binding_kcal),
            func.sum(Result.compute_seconds),
        ).where(Result.run_id == run_id)
    ).one()

    return {
        "run_id": run_id,
        "n_results": row[0],
        "n_variants": row[1],
        "n_conditions": row[2],
        "delta_g_binding_min": row[3],
        "delta_g_binding_max": row[4],
        "compute_seconds": row[5],
    }


def get_or_create_target(
    session: Session, name: str, sequence: str, description: str | None = None
) -> Target:
    """Fetch the named target, creating it if this is the first sweep against it."""
    target = session.execute(select(Target).where(Target.name == name)).scalar_one_or_none()
    if target is None:
        target = Target(name=name, sequence=sequence.strip().upper(), description=description)
        session.add(target)
        session.flush()
    elif target.sequence != sequence.strip().upper():
        raise ValueError(
            f"target {name} is already in the database with a different sequence. "
            "Give the new target its own name so existing results stay meaningful."
        )
    return target
