"""Catalog loading: derived columns, idempotency, and the generated DDL check."""

from __future__ import annotations

import subprocess
import sys

import pytest
from sqlalchemy import func, select

from nupack_campaign.catalog import CatalogError, classify, load_variants, read_catalog
from nupack_campaign.models import Variant

HEADER = (
    "sequence,length,gc_content,trunc_5prime,trunc_3prime,"
    "position_1,original_base_1,mutated_base_1,mutation_type_1,"
    "position_2,original_base_2,mutated_base_2,mutation_type_2,"
    "num_transitions,num_transversions\n"
)


def write_csv(tmp_path, rows: str, header: str = HEADER):
    path = tmp_path / "scan.csv"
    path.write_text(header + rows)
    return path


def test_generated_sql_matches_the_models(repo_root):
    result = subprocess.run(
        [sys.executable, "scripts/dump_schema_sql.py", "--check"],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_classify_covers_both_axes():
    assert classify(0, 0) == "reference"
    assert classify(0, 3) == "truncation_only"
    assert classify(1, 0) == "single_mutation"
    assert classify(2, 0) == "double_mutation"
    assert classify(2, 4) == "double_mutation_truncated"


def test_derived_columns_come_from_the_sequence(tmp_path):
    path = write_csv(tmp_path, "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n")
    catalog = read_catalog(path)
    row = catalog.iloc[0]

    assert row["length_nt"] == 8
    assert row["gc_content"] == pytest.approx(0.5)
    assert row["num_mutations"] == 0
    assert row["variant_class"] == "reference"
    assert row["name"] == "t0_0_wt"


def test_mutation_signature_and_class(tmp_path):
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,1,2,0,A,G,transition,3,T,C,transition,2,0\n",
    )
    row = read_catalog(path).iloc[0]

    assert row["num_mutations"] == 2
    assert row["trunc_total"] == 3
    assert row["mutation_signature"] == "A0G,T3C"
    assert row["variant_class"] == "double_mutation_truncated"
    assert row["name"] == "t1_2_A0G,T3C"


def test_a_bare_sequence_column_is_enough(tmp_path):
    path = write_csv(tmp_path, "GGCCAATT\nAATTGGCC\n", header="sequence\n")
    catalog = read_catalog(path)
    assert len(catalog) == 2
    assert set(catalog["variant_class"]) == {"reference"}


def test_length_column_is_cross_checked(tmp_path):
    path = write_csv(tmp_path, "GGCCAATT,99,0.5,0,0,,,,none,,,,none,0,0\n")
    with pytest.raises(CatalogError, match="length column says 99"):
        read_catalog(path)


def test_gc_column_is_cross_checked(tmp_path):
    path = write_csv(tmp_path, "GGCCAATT,8,0.99,0,0,,,,none,,,,none,0,0\n")
    with pytest.raises(CatalogError, match="gc_content column says"):
        read_catalog(path)


def test_degenerate_bases_are_rejected(tmp_path):
    path = write_csv(tmp_path, "GGCCAANN,8,0.5,0,0,,,,none,,,,none,0,0\n")
    with pytest.raises(CatalogError, match="non-ACGT"):
        read_catalog(path)


def test_colliding_scan_coordinates_fall_back_to_a_sequence_hash(tmp_path):
    """Two different sequences at the same coordinates still get unique names."""
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"
        "AATTGGCC,8,0.5,0,0,,,,none,,,,none,0,0\n",
    )
    names = list(read_catalog(path)["name"])

    assert len(set(names)) == 2
    assert all(name.startswith("t0_0_wt_") for name in names)


def test_the_hash_suffix_is_stable_across_reloads(tmp_path):
    rows = "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\nAATTGGCC,8,0.5,0,0,,,,none,,,,none,0,0\n"
    first = list(read_catalog(write_csv(tmp_path, rows))["name"])

    reordered = tmp_path / "reordered.csv"
    reordered.write_text(HEADER + "".join(reversed(rows.splitlines(keepends=True))))
    second = list(read_catalog(reordered)["name"])

    assert sorted(first) == sorted(second)


def test_a_genuinely_duplicated_row_is_rejected(tmp_path):
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n",
    )
    with pytest.raises(CatalogError, match="does not uniquely identify"):
        read_catalog(path)


def test_an_explicit_name_column_wins(tmp_path):
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n",
        header="name," + HEADER,
    )
    path.write_text(
        "name,sequence,length,gc_content,trunc_5prime,trunc_3prime,position_1,"
        "original_base_1,mutated_base_1,mutation_type_1,position_2,original_base_2,"
        "mutated_base_2,mutation_type_2,num_transitions,num_transversions\n"
        "probe_A,GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"
    )
    assert read_catalog(path)["name"].iloc[0] == "probe_A"


def test_loading_is_idempotent(session, tmp_path):
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"
        "GGCCAATT,8,0.5,1,0,,,,none,,,,none,0,0\n",
    )
    first = load_variants(session, path)
    second = load_variants(session, path)

    assert (first.inserted, second.inserted) == (2, 0)
    assert second.skipped == 2
    assert session.execute(select(func.count(Variant.variant_id))).scalar_one() == 2


def test_a_changed_sequence_under_an_existing_name_is_rejected(session, tmp_path):
    load_variants(session, write_csv(tmp_path, "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"))

    clash = tmp_path / "other.csv"
    clash.write_text(HEADER + "AATTGGCC,8,0.5,0,0,,,,none,,,,none,0,0\n")
    with pytest.raises(CatalogError, match="different"):
        load_variants(session, clash)


def test_load_stats_report_the_class_breakdown(session, tmp_path):
    path = write_csv(
        tmp_path,
        "GGCCAATT,8,0.5,0,0,,,,none,,,,none,0,0\n"
        "GCCAATT,7,0.4286,1,0,,,,none,,,,none,0,0\n"
        "GGCCAAGG,8,0.75,0,0,6,T,G,transversion,7,T,G,transversion,0,2\n",
    )
    stats = load_variants(session, path)
    assert stats.by_class == {"reference": 1, "truncation_only": 1, "double_mutation": 1}
