"""Ranking and trend behaviour.

The direction tests exist because the sort order is the whole point of the
query and a well meaning ASC would look like a tidy-up rather than a change in
meaning.
"""

from __future__ import annotations

import dataclasses

import pytest
from sqlalchemy import func, select

from nupack_campaign.models import Selection, SelectionMember, Variant
from nupack_campaign.selection import (
    ConditionNotFound,
    SelectionCriteria,
    explain_query_plan,
    persist_selection,
    resolve_condition_id,
    select_candidates,
)
from nupack_campaign.sweep import expand_grid, get_or_create_target, run_sweep, start_run
from nupack_campaign.trends import DIMENSIONS, UnknownDimension, resolve, trend_by

from .conftest import fixture_name

WIDE_OPEN = SelectionCriteria(
    temperature_c=37.0,
    na_molar=0.15,
    mg_molar=0.0,
    material="dna04",
    ensemble="stacking",
    min_delta_g_binding=-999.0,
    max_delta_g_binding=999.0,
    min_ddg=-999.0,
    gc_min=0.0,
    gc_max=1.0,
    trunc_total_min=0,
    trunc_total_max=99,
    num_mutations_min=0,
    num_mutations_max=2,
    per_class_cap=99,
    n_final=99,
)


@pytest.fixture
def swept(populated_session, thermo_engine, grid, target):
    run = start_run(populated_session, thermo_engine, grid)
    populated_session.commit()
    run_sweep(
        populated_session, thermo_engine, run, target, expand_grid(grid), verbose=False
    )
    return populated_session, run


def rank(session, run, criteria, query_path):
    return select_candidates(session, run.run_id, criteria, query_path=query_path)


# --- ranking direction -------------------------------------------------------


def test_ranking_is_most_positive_first(swept, query_path):
    """The headline behaviour. Do not let this flip to ASC."""
    session, run = swept
    shortlist = rank(session, run, WIDE_OPEN, query_path)

    assert len(shortlist) == 9
    assert shortlist["delta_g_binding_kcal"].is_monotonic_decreasing
    assert shortlist["delta_g_binding_kcal"].iloc[0] == shortlist["delta_g_binding_kcal"].max()


def test_strongest_binders_query_is_the_mirror(swept, query_path, strongest_query_path):
    session, run = swept
    disrupted = rank(session, run, WIDE_OPEN, query_path)
    strongest = rank(session, run, WIDE_OPEN, strongest_query_path)

    assert strongest["delta_g_binding_kcal"].is_monotonic_increasing
    assert list(strongest["variant_id"]) == list(reversed(disrupted["variant_id"]))


def test_the_two_queries_take_the_same_parameters(query_path, strongest_query_path):
    """They are two files only because ORDER BY cannot take a bind parameter."""
    import re

    def markers(path):
        body = "\n".join(
            line for line in path.read_text().splitlines() if not line.strip().startswith("--")
        )
        return set(re.findall(r":(\w+)", body))

    assert markers(query_path) == markers(strongest_query_path)


def test_top_of_the_ranking_is_a_mutant_not_the_reference(swept, query_path):
    """Mutations cost binding, so the untruncated reference should sink."""
    session, run = swept
    shortlist = rank(session, run, WIDE_OPEN, query_path)
    assert shortlist["num_mutations"].iloc[0] > 0
    assert shortlist.iloc[-1]["name"] == fixture_name(0, 0, 0)


# --- ddG against the per-truncation reference --------------------------------


def test_reference_rows_have_zero_ddg(swept, query_path):
    session, run = swept
    shortlist = rank(session, run, WIDE_OPEN, query_path)
    unmutated = shortlist[shortlist.num_mutations == 0]

    assert len(unmutated) == 3
    assert (unmutated["ddg_vs_reference_kcal"].abs() < 1e-9).all()


def test_ddg_uses_the_matching_truncation_not_the_global_wild_type(swept, query_path):
    """Each truncation gets its own baseline, so truncation cost is not in ddG."""
    session, run = swept
    shortlist = rank(session, run, WIDE_OPEN, query_path).set_index("name")

    for trunc_5, trunc_3 in ((0, 0), (1, 0), (1, 1)):
        reference = shortlist.loc[
            fixture_name(trunc_5, trunc_3, 0), "delta_g_binding_kcal"
        ]
        mutant = shortlist.loc[fixture_name(trunc_5, trunc_3, 1)]

        assert mutant["reference_delta_g_binding_kcal"] == pytest.approx(reference)
        # The two columns are rounded independently in the query, so compare to
        # the precision they are stored at rather than to full float precision.
        assert mutant["ddg_vs_reference_kcal"] == pytest.approx(
            mutant["delta_g_binding_kcal"] - reference, abs=0.002
        )


def test_min_ddg_filters_out_the_references(swept, query_path):
    session, run = swept
    criteria = dataclasses.replace(WIDE_OPEN, min_ddg=0.001)
    shortlist = rank(session, run, criteria, query_path)
    assert (shortlist["num_mutations"] > 0).all()


# --- filters and caps --------------------------------------------------------


def test_per_class_cap_spreads_the_shortlist(swept, query_path):
    session, run = swept
    shortlist = rank(session, run, dataclasses.replace(WIDE_OPEN, per_class_cap=1), query_path)
    assert len(shortlist) == shortlist["variant_class"].nunique()


def test_n_final_caps_the_shortlist(swept, query_path):
    session, run = swept
    assert len(rank(session, run, dataclasses.replace(WIDE_OPEN, n_final=4), query_path)) == 4


def test_truncation_filter_is_applied(swept, query_path):
    session, run = swept
    shortlist = rank(
        session, run, dataclasses.replace(
            WIDE_OPEN,
            trunc_total_min=1,
            trunc_total_max=1,
        ), query_path
    )
    assert (shortlist["trunc_total"] == 1).all()


def test_mutation_count_filter_is_applied(swept, query_path):
    session, run = swept
    shortlist = rank(
        session, run, dataclasses.replace(
            WIDE_OPEN,
            num_mutations_min=0,
            num_mutations_max=0,
        ), query_path
    )
    assert (shortlist["num_mutations"] == 0).all()


def test_gc_filter_is_applied(swept, query_path):
    session, run = swept
    shortlist = rank(
        session,
        run,
        dataclasses.replace(WIDE_OPEN,
        gc_min=0.6,
        gc_max=0.7),
        query_path,
    )
    assert shortlist["gc_content"].between(0.6, 0.7).all()


def test_binding_window_is_applied(swept, query_path):
    session, run = swept
    wide = rank(session, run, WIDE_OPEN, query_path)
    centre = float(wide["delta_g_binding_kcal"].sort_values().iloc[len(wide) // 2])
    narrow = rank(
        session,
        run,
        dataclasses.replace(
            WIDE_OPEN,
            min_delta_g_binding=centre - 0.5,
            max_delta_g_binding=centre + 0.5,
        ),
        query_path,
    )
    assert 0 < len(narrow) < len(wide)


def test_empty_shortlist_still_has_its_columns(swept, query_path):
    session, run = swept
    criteria = dataclasses.replace(WIDE_OPEN, min_delta_g_binding=500.0, max_delta_g_binding=600.0)
    shortlist = rank(session, run, criteria, query_path)
    assert shortlist.empty
    assert "delta_g_binding_kcal" in shortlist.columns
    assert "ddg_vs_reference_kcal" in shortlist.columns


def test_unknown_condition_raises(swept, query_path):
    session, run = swept
    with pytest.raises(ConditionNotFound, match="99.0 C"):
        rank(session, run, dataclasses.replace(WIDE_OPEN, temperature_c=99.0), query_path)


def test_condition_resolution_tolerates_float_representation(swept):
    session, _ = swept
    criteria = dataclasses.replace(WIDE_OPEN, na_molar=0.1 + 0.05)
    assert resolve_condition_id(session, criteria) > 0


# --- persistence -------------------------------------------------------------


def test_persist_selection_stores_members_and_ddg(swept, query_path):
    session, run = swept
    shortlist = rank(session, run, WIDE_OPEN, query_path)
    record = persist_selection(session, run.run_id, WIDE_OPEN, shortlist, query_path=query_path)

    assert record.n_selected == len(shortlist)
    assert record.n_candidates_in == 9

    members = session.execute(
        select(SelectionMember).where(SelectionMember.selection_id == record.selection_id)
    ).scalars().all()
    assert sorted(m.rank for m in members) == list(range(1, len(shortlist) + 1))
    assert all(m.ddg_vs_reference_kcal is not None for m in members)


def test_two_selections_on_one_run_are_both_kept(swept, query_path):
    session, run = swept
    for cap in (1, 99):
        criteria = dataclasses.replace(WIDE_OPEN, per_class_cap=cap)
        persist_selection(
            session, run.run_id, criteria, rank(
                session,
                run,
                criteria,
                query_path,
            ), query_path=query_path
        )
    assert session.execute(select(func.count(Selection.selection_id))).scalar_one() == 2


# --- trends ------------------------------------------------------------------


def test_every_dimension_runs(swept, trend_query_path):
    session, run = swept
    for name in DIMENSIONS:
        frame = trend_by(session, run.run_id, name, query_path=trend_query_path)
        assert "n" in frame.columns
        assert "mean_delta_g_binding" in frame.columns


def test_trend_by_truncation_has_one_row_per_level(swept, trend_query_path):
    session, run = swept
    frame = trend_by(session, run.run_id, "truncation", query_path=trend_query_path)
    assert list(frame["bases trimmed"]) == [0, 1, 2]
    assert frame["n"].sum() == 9 * 2  # nine variants at two conditions


def test_trend_group_counts_match_the_catalog(swept, trend_query_path):
    session, run = swept
    frame = trend_by(session, run.run_id, "num_mutations", query_path=trend_query_path)
    counts = dict(zip(frame["mutations"], frame["n"], strict=True))

    expected = dict(
        session.execute(
            select(Variant.num_mutations, func.count(Variant.variant_id)).group_by(
                Variant.num_mutations
            )
        ).all()
    )
    assert counts == {key: value * 2 for key, value in expected.items()}


def test_trend_can_be_scoped_to_one_condition(swept, trend_query_path):
    session, run = swept
    condition_id = resolve_condition_id(session, WIDE_OPEN)
    scoped = trend_by(
        session,
        run.run_id,
        "truncation",
        condition_id=condition_id,
        query_path=trend_query_path,
    )
    pooled = trend_by(session, run.run_id, "truncation", query_path=trend_query_path)
    assert scoped["n"].sum() * 2 == pooled["n"].sum()


def test_min_group_size_drops_small_buckets(swept, trend_query_path):
    session, run = swept
    frame = trend_by(
        session,
        run.run_id,
        "truncation",
        min_group_size=1000,
        query_path=trend_query_path,
    )
    assert frame.empty


def test_trend_reference_rows_have_zero_mean_ddg(swept, trend_query_path):
    session, run = swept
    frame = trend_by(session, run.run_id, "num_mutations", query_path=trend_query_path)
    unmutated = frame[frame["mutations"] == 0]
    assert unmutated["mean_ddg_vs_reference"].iloc[0] == pytest.approx(0.0, abs=1e-6)


def test_unknown_dimension_is_rejected():
    with pytest.raises(UnknownDimension, match="not a known dimension"):
        resolve("'; DROP TABLE results; --")


def test_dimension_expressions_only_touch_known_aliases():
    """The expressions are interpolated into SQL, so keep them to the two
    aliases the template defines and free of anything that could close the
    statement."""
    import re

    for name, dimension in DIMENSIONS.items():
        aliases = set(re.findall(r"\b([a-z]+)\.", dimension.expression))
        assert aliases <= {"v", "c", "a"}, f"{name} references an unknown alias"
        assert ";" not in dimension.expression
        assert "--" not in dimension.expression


# --- query plan --------------------------------------------------------------


def test_query_plan_uses_the_ranking_index(swept, query_path):
    session, run = swept
    plan = explain_query_plan(session, run.run_id, WIDE_OPEN, query_path=query_path)
    assert "idx_results_ranking" in plan


def test_repo_selection_config_is_valid(repo_root):
    SelectionCriteria.from_yaml(repo_root / "config" / "selection.yaml")


def test_criteria_validation_catches_an_inverted_band():
    with pytest.raises(ValueError, match="band is empty"):
        dataclasses.replace(WIDE_OPEN, min_delta_g_binding=5.0, max_delta_g_binding=-5.0).validate()


def test_get_or_create_target_rejects_a_changed_sequence(populated_session):
    get_or_create_target(populated_session, "t", "ACGT")
    populated_session.commit()
    with pytest.raises(ValueError, match="different sequence"):
        get_or_create_target(populated_session, "t", "TGCA")
