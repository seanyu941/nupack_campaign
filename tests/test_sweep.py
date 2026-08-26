"""Sweep behaviour: grid expansion, condition reuse, resumability, and the
three-body binding identity."""

from __future__ import annotations

import pytest
from sqlalchemy import func, select

from nupack_campaign.engines import ConditionPoint, StubEngine
from nupack_campaign.models import Condition, Result, TargetEnergy
from nupack_campaign.sweep import (
    ensure_conditions,
    expand_grid,
    grid_hash,
    load_grid_file,
    pending_pairs,
    resume_run,
    run_sweep,
    start_run,
)


def sweep_once(session, thermo_engine, run, target, grid):
    return run_sweep(session, thermo_engine, run, target, expand_grid(grid), verbose=False)


def test_expand_grid_is_the_cross_product(grid):
    points = expand_grid(grid)
    assert len(points) == 2
    assert len({point.as_key() for point in points}) == len(points)


def test_grid_hash_ignores_key_order(grid):
    assert grid_hash(grid) == grid_hash(dict(reversed(list(grid.items()))))


def test_grid_hash_changes_when_a_value_is_added(grid):
    widened = {**grid, "temperature_c": grid["temperature_c"] + [45.0]}
    assert grid_hash(widened) != grid_hash(grid)


def test_repo_grid_file_loads(repo_root):
    assert expand_grid(load_grid_file(repo_root / "config" / "sweep_grid.yaml"))


def test_ensure_conditions_reuses_existing_rows(session, grid):
    points = expand_grid(grid)
    first = ensure_conditions(session, points)
    session.commit()
    second = ensure_conditions(session, points)
    session.commit()

    assert first == second
    assert session.execute(select(func.count(Condition.condition_id))).scalar_one() == len(points)


def test_sweep_fills_every_cell(populated_session, thermo_engine, grid, target):
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    stats = sweep_once(populated_session, thermo_engine, run, target, grid)

    assert stats.total_cells == 9 * 2
    assert stats.computed_cells == stats.total_cells
    assert run.status == "complete"


def test_rerunning_a_finished_sweep_computes_nothing(
    populated_session,
    thermo_engine,
    grid,
    target,
):
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    sweep_once(populated_session, thermo_engine, run, target, grid)

    before = populated_session.execute(select(func.count(Result.result_id))).scalar_one()
    second = sweep_once(populated_session, thermo_engine, run, target, grid)
    after = populated_session.execute(select(func.count(Result.result_id))).scalar_one()

    assert second.computed_cells == 0
    assert before == after


def test_widening_the_grid_only_computes_the_difference(
    populated_session,
    thermo_engine,
    grid,
    target,
):
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    sweep_once(populated_session, thermo_engine, run, target, grid)

    widened = {**grid, "temperature_c": grid["temperature_c"] + [45.0]}
    second = sweep_once(populated_session, thermo_engine, run, target, widened)
    assert second.computed_cells == 9


def test_target_energy_is_computed_once_per_condition(
    populated_session,
    thermo_engine,
    grid,
    target,
):
    """The whole point of the target_energies table."""
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    sweep_once(populated_session, thermo_engine, run, target, grid)

    n_target_rows = populated_session.execute(
        select(func.count(TargetEnergy.target_energy_id))
    ).scalar_one()
    n_results = populated_session.execute(select(func.count(Result.result_id))).scalar_one()

    assert n_target_rows == 2       # one per condition
    assert n_results == 18          # nine variants at two conditions


def test_binding_identity_holds_for_every_row(populated_session, thermo_engine, grid, target):
    """delta_g_binding must equal complex - cdna - target, or the column lies."""
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    sweep_once(populated_session, thermo_engine, run, target, grid)

    target_dg = {
        row.condition_id: row.free_energy_kcal
        for row in populated_session.execute(select(TargetEnergy)).scalars()
    }

    for result in populated_session.execute(select(Result)).scalars():
        implied = (
            result.complex_dg_kcal - result.cdna_dg_kcal - target_dg[result.condition_id]
        )
        assert result.delta_g_binding_kcal == pytest.approx(implied, abs=1e-3)


def test_pending_pairs_empties_out_as_the_sweep_runs(
    populated_session,
    thermo_engine,
    grid,
    target,
):
    points = expand_grid(grid)
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()

    ensure_conditions(populated_session, points)
    populated_session.commit()
    assert len(pending_pairs(populated_session, run.run_id, target.target_id)) == 18

    sweep_once(populated_session, thermo_engine, run, target, grid)
    assert pending_pairs(populated_session, run.run_id, target.target_id) == []


def test_resume_rejects_a_changed_grid(populated_session, thermo_engine, grid):
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    with pytest.raises(ValueError, match="different sweep grid"):
        resume_run(populated_session, run.run_id, {**grid, "na_molar": [0.15, 1.0]})


def test_stub_engine_is_deterministic():
    point = ConditionPoint(37.0, 0.15, 0.0, "dna04", "stacking")
    first = StubEngine().evaluate_binding("ACGTACGTACGT", "ACGTACGTACGT", point, -12.0)
    second = StubEngine().evaluate_binding("ACGTACGTACGT", "ACGTACGTACGT", point, -12.0)
    assert first == second


def test_ensemble_free_energy_never_exceeds_mfe():
    engine = StubEngine()
    point = ConditionPoint(37.0, 0.15, 0.0, "dna04", "stacking")
    for sequence in ("ACGTACGTACGT", "GGCCGGCCGGCC", "ATATATATATAT"):
        result = engine.evaluate_binding(sequence, "ACGTACGTACGT", point, -12.0)
        assert result.cdna_dg_kcal <= result.mfe_kcal


def test_a_matching_design_binds_better_than_a_mismatched_one():
    """Sanity check on the stub, so the ranking tests are testing something."""
    engine = StubEngine()
    point = ConditionPoint(37.0, 0.15, 0.0, "dna04", "stacking")
    target_sequence = "AAGCGCGGAAGC"
    perfect = "GCTTCCGCGCTT"
    mismatched = "GAAACCGCGCAA"

    good = engine.evaluate_binding(perfect, target_sequence, point, -12.0)
    bad = engine.evaluate_binding(mismatched, target_sequence, point, -12.0)
    assert good.delta_g_binding_kcal < bad.delta_g_binding_kcal
