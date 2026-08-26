"""Target alignment.

A design is much shorter than its target, so where it binds has to be found
rather than assumed. These tests pin the geometry, because getting it wrong
produces results that look completely normal and mean nothing.
"""

from __future__ import annotations

import pytest
from sqlalchemy import func, select

from nupack_campaign.catalog import align_variants
from nupack_campaign.engines import (
    NoBindingSite,
    TargetAligner,
    check_binding_site,
    reverse_complement,
)
from nupack_campaign.models import VariantAlignment

from .conftest import BINDING_OFFSET, PARENT, TARGET_SEQUENCE, fixture_name


def test_reverse_complement_round_trips():
    assert reverse_complement("GCTTCC") == "GGAAGC"
    assert reverse_complement(reverse_complement(PARENT)) == PARENT


def test_a_perfect_design_finds_its_window():
    alignment = TargetAligner(TARGET_SEQUENCE).align(PARENT)
    assert alignment.offset == BINDING_OFFSET
    assert alignment.n_mismatch == 0
    assert alignment.match_fraction == 1.0


def test_the_window_is_not_assumed_to_start_at_zero():
    """The whole reason the aligner exists."""
    assert TARGET_SEQUENCE.find(reverse_complement(PARENT)) == BINDING_OFFSET
    assert BINDING_OFFSET > 0


def test_each_mutation_costs_one_paired_base():
    aligner = TargetAligner(TARGET_SEQUENCE)
    perfect = aligner.align(PARENT)

    mutated = list(PARENT)
    mutated[2] = "A" if mutated[2] != "A" else "C"
    one_off = aligner.align("".join(mutated))

    assert one_off.n_complementary == perfect.n_complementary - 1


def test_truncating_the_3_prime_end_moves_the_window():
    """Design and target run antiparallel, so trimming the 3' end shifts the
    window along the target rather than shortening it in place."""
    aligner = TargetAligner(TARGET_SEQUENCE)
    full = aligner.align(PARENT)
    trimmed = aligner.align(PARENT[:-2])

    assert trimmed.offset == full.offset + 2
    assert trimmed.n_mismatch == 0


def test_truncating_the_5_prime_end_leaves_the_offset_alone():
    aligner = TargetAligner(TARGET_SEQUENCE)
    full = aligner.align(PARENT)
    trimmed = aligner.align(PARENT[2:])

    assert trimmed.offset == full.offset
    assert trimmed.length == full.length - 2
    assert trimmed.n_mismatch == 0


def test_alignment_is_cached_but_not_confused():
    aligner = TargetAligner(TARGET_SEQUENCE)
    first = aligner.align(PARENT)
    other = aligner.align(PARENT[2:])
    again = aligner.align(PARENT)

    assert first == again
    assert other != first


def test_a_design_longer_than_its_target_still_aligns():
    aligner = TargetAligner("ACGT")
    alignment = aligner.align("ACGTACGTACGT")
    assert alignment.offset == 0


# --- the wrong target guard ---------------------------------------------------


def test_check_binding_site_accepts_the_real_target():
    site = check_binding_site([PARENT], TARGET_SEQUENCE)
    assert site.offset == BINDING_OFFSET


def test_check_binding_site_rejects_an_unrelated_target():
    with pytest.raises(NoBindingSite, match="right strand"):
        check_binding_site([PARENT], "TTTTTTTTTTTTTTTTTTTTTTTT")


def test_check_binding_site_rejects_the_target_given_the_wrong_way_round():
    """Passing the design's own strand instead of what it binds is the easy
    mistake, and it is not caught by anything else in the pipeline."""
    with pytest.raises(NoBindingSite):
        check_binding_site([PARENT], PARENT + PARENT)


def test_check_binding_site_looks_at_the_best_design_not_every_design():
    """A mutation scan is meant to contain designs that bind badly, so a single
    poor design must not fail the whole import."""
    junk = "TTTTTTTTTTTT"
    assert check_binding_site([junk, PARENT], TARGET_SEQUENCE).n_mismatch == 0


# --- the stored table ---------------------------------------------------------


def test_align_variants_writes_one_row_per_design(populated_session, target):
    written = align_variants(populated_session, target)
    total = populated_session.execute(
        select(func.count(VariantAlignment.variant_id))
    ).scalar_one()

    assert written == 9
    assert total == 9


def test_align_variants_is_idempotent(populated_session, target):
    align_variants(populated_session, target)
    second = align_variants(populated_session, target)

    assert second == 0
    total = populated_session.execute(
        select(func.count(VariantAlignment.variant_id))
    ).scalar_one()
    assert total == 9


def test_stored_alignment_matches_the_aligner(populated_session, aligned_target):
    from nupack_campaign.models import Variant

    aligner = TargetAligner(aligned_target.sequence)
    rows = populated_session.execute(
        select(Variant.sequence, VariantAlignment).join(
            VariantAlignment, VariantAlignment.variant_id == Variant.variant_id
        )
    ).all()

    assert len(rows) == 9
    for sequence, stored in rows:
        expected = aligner.align(sequence)
        assert stored.target_offset == expected.offset
        assert stored.n_complementary == expected.n_complementary
        assert stored.n_mismatch == expected.n_mismatch


def test_paired_bases_folds_truncation_and_mutation_together(populated_session, aligned_target):
    """The point of storing n_complementary as its own column."""
    from nupack_campaign.models import Variant

    by_name = {
        name: alignment
        for name, alignment in populated_session.execute(
            select(Variant.name, VariantAlignment).join(
                VariantAlignment, VariantAlignment.variant_id == Variant.variant_id
            )
        ).all()
    }

    # Untruncated and unmutated pairs everywhere.
    assert by_name[fixture_name(0, 0, 0)].n_complementary == len(PARENT)
    # One base trimmed costs one pair, two mutations cost two.
    assert by_name[fixture_name(1, 0, 0)].n_complementary == len(PARENT) - 1
    assert by_name[fixture_name(0, 0, 1)].n_complementary == len(PARENT) - 2
    # Both axes together: two trimmed plus two mutated.
    assert by_name[fixture_name(1, 1, 1)].n_complementary == len(PARENT) - 4
