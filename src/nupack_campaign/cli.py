"""Command line entry point.

    python -m nupack_campaign.cli load-variants --catalog scan.csv
    python -m nupack_campaign.cli import-results --catalog scan.csv --target-sequence ...
    python -m nupack_campaign.cli sweep --engine stub --target-sequence ...
    python -m nupack_campaign.cli rank --run-id 1
    python -m nupack_campaign.cli trends --by truncation
    python -m nupack_campaign.cli explain --run-id 1

The Makefile wraps these. ``import-results`` is the quick route on a scan CSV
that already has free energies in it: it unpivots the wide energy columns into
``results`` so the ranking and trend queries work without rerunning NUPACK.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from sqlalchemy import desc, select

from .catalog import load_variants
from .db import DEFAULT_DB_URL, init_db, make_engine, session_scope
from .engines import get_engine
from .importer import import_scan_csv
from .models import Condition, Run, Selection
from .selection import (
    STRONGEST_QUERY_PATH,
    SelectionCriteria,
    explain_query_plan,
    persist_selection,
    resolve_condition_id,
    select_candidates,
)
from .sweep import (
    expand_grid,
    get_or_create_target,
    load_grid_file,
    resume_run,
    run_summary,
    run_sweep,
    start_run,
)
from .trends import DIMENSIONS, describe_dimensions, trend_by

DEFAULT_GRID = "config/sweep_grid.yaml"
DEFAULT_SELECTION = "config/selection.yaml"
DEFAULT_CATALOG = "data/variants.csv"
DEFAULT_QUERY = "sql/select_candidates.sql"

RANK_COLUMNS = [
    "name",
    "variant_class",
    "mutation_signature",
    "trunc_5prime",
    "trunc_3prime",
    "length_nt",
    "gc_content",
    "delta_g_binding_kcal",
    "ddg_vs_reference_kcal",
    "rank_in_class",
]


def latest_complete_run_id(session) -> int:
    run_id = session.execute(
        select(Run.run_id).where(Run.status == "complete").order_by(desc(Run.run_id)).limit(1)
    ).scalar_one_or_none()
    if run_id is None:
        raise SystemExit("No completed run found. Run the sweep or import first.")
    return int(run_id)


def show(frame: pd.DataFrame, columns: list[str] | None = None) -> None:
    if frame.empty:
        print("(no rows)")
        return
    if columns:
        frame = frame[[c for c in columns if c in frame.columns]]
    with pd.option_context("display.width", 170, "display.max_columns", 30):
        print(frame.to_string(index=False))


def cmd_init_db(args: argparse.Namespace) -> int:
    init_db(make_engine(args.db))
    print(f"[init-db] schema ready at {args.db}")
    return 0


def cmd_load_variants(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    init_db(engine)
    with session_scope(engine) as session:
        stats = load_variants(session, args.catalog)
    print(f"[load-variants] {stats}")
    return 0


def cmd_import_results(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    init_db(engine)
    with session_scope(engine) as session:
        if args.load_variants:
            stats = load_variants(session, args.catalog)
            print(f"[load-variants] {stats}")

        result = import_scan_csv(
            session,
            args.catalog,
            target_name=args.target_name,
            target_sequence=args.target_sequence,
            na_molar=args.na_molar,
            mg_molar=args.mg_molar,
            note=args.note,
        )

    print(f"[import] {result}")
    print(
        f"[import] best binding window starts at target offset "
        f"{result.binding_offset}, {result.binding_paired} bases paired"
    )
    for temperature, energy in sorted(result.target_energies.items()):
        print(f"[import] recovered target dG at {temperature:g}C: {energy:.4f} kcal/mol")
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    init_db(engine)

    grid = load_grid_file(args.grid)
    points = expand_grid(grid)
    thermo = get_engine(args.engine)

    with session_scope(engine) as session:
        target = get_or_create_target(session, args.target_name, args.target_sequence)
        run = (
            resume_run(session, args.resume, grid)
            if args.resume
            else start_run(session, thermo, grid, note=args.note)
        )
        session.commit()

        stats = run_sweep(
            session,
            thermo,
            run,
            target,
            points,
            batch_size=args.batch_size,
            verbose=not args.quiet,
        )
        summary = run_summary(session, stats.run_id)

    print(
        f"[sweep] run {stats.run_id}: {summary['n_results']} rows, "
        f"delta_g_binding from {summary['delta_g_binding_min']:.2f} to "
        f"{summary['delta_g_binding_max']:.2f} kcal/mol"
    )
    print(f"[sweep] computed {stats.computed_cells}, reused {stats.skipped_cells}")
    return 0


def cmd_rank(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    criteria = SelectionCriteria.from_yaml(args.config)
    query_path = Path(STRONGEST_QUERY_PATH) if args.strongest else Path(args.query)

    with session_scope(engine) as session:
        run_id = args.run_id or latest_complete_run_id(session)
        shortlist = select_candidates(session, run_id, criteria, query_path=query_path)

        if shortlist.empty:
            print("[rank] nothing passed the filters, loosen config/selection.yaml")
            return 1

        record = persist_selection(session, run_id, criteria, shortlist, query_path=query_path)
        direction = "most negative first" if args.strongest else "most positive first"
        print(
            f"[rank] selection {record.selection_id}: {record.n_candidates_in} candidates -> "
            f"{record.n_selected} selected, ranked by delta_g_binding {direction}"
        )

    show(shortlist, RANK_COLUMNS)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        shortlist.to_csv(args.out, index=False)
        print(f"[rank] written to {args.out}")
    return 0


def cmd_trends(args: argparse.Namespace) -> int:
    if args.list:
        show(describe_dimensions())
        return 0

    engine = make_engine(args.db)
    with session_scope(engine) as session:
        run_id = args.run_id or latest_complete_run_id(session)

        condition_id = None
        if not args.all_conditions:
            criteria = SelectionCriteria.from_yaml(args.config)
            condition_id = resolve_condition_id(session, criteria)

        for dimension in args.by:
            frame = trend_by(
                session,
                run_id,
                dimension,
                condition_id=condition_id,
                min_group_size=args.min_group_size,
            )
            print(f"\n[trends] {dimension} ({DIMENSIONS[dimension].description})")
            show(frame)

            if args.out:
                out = Path(args.out)
                out.mkdir(parents=True, exist_ok=True)
                frame.to_csv(out / f"trend_{dimension}.csv", index=False)

    if args.out:
        print(f"\n[trends] written to {args.out}/")
    return 0


def cmd_explain(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    criteria = SelectionCriteria.from_yaml(args.config)
    with session_scope(engine) as session:
        run_id = args.run_id or latest_complete_run_id(session)
        print(explain_query_plan(session, run_id, criteria, query_path=Path(args.query)))
    return 0


def cmd_runs(args: argparse.Namespace) -> int:
    engine = make_engine(args.db)
    with session_scope(engine) as session:
        for run in session.execute(select(Run).order_by(Run.run_id)).scalars():
            summary = run_summary(session, run.run_id)
            print(
                f"run {run.run_id:>3}  {run.status:<9} engine={run.engine_name:<9} "
                f"rows={summary['n_results']:<7} conditions={summary['n_conditions']:<3} "
                f"git={run.git_sha[:8]}"
            )
        conditions = session.execute(
            select(Condition).order_by(Condition.condition_id)
        ).scalars()
        for condition in conditions:
            print(
                f"  condition {condition.condition_id:>3}: {condition.temperature_c:g}C "
                f"Na={condition.na_molar:g}M Mg={condition.mg_molar:g}M"
            )
        for record in session.execute(select(Selection).order_by(Selection.selection_id)).scalars():
            print(
                f"selection {record.selection_id:>3} on run {record.run_id}: "
                f"{record.n_candidates_in} -> {record.n_selected}"
            )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nupack_campaign", description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB_URL, help="SQLAlchemy database URL")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init-db", help="create tables and indexes")
    init.set_defaults(func=cmd_init_db)

    load = subparsers.add_parser("load-variants", help="load a design scan CSV")
    load.add_argument("--catalog", default=DEFAULT_CATALOG)
    load.set_defaults(func=cmd_load_variants)

    imp = subparsers.add_parser(
        "import-results", help="unpivot a scan CSV that already has free energies"
    )
    imp.add_argument("--catalog", default=DEFAULT_CATALOG)
    imp.add_argument("--target-name", default="target")
    imp.add_argument("--target-sequence", required=True)
    imp.add_argument("--na-molar", type=float, default=0.15)
    imp.add_argument("--mg-molar", type=float, default=0.0)
    imp.add_argument("--note", default=None)
    imp.add_argument(
        "--no-load-variants",
        dest="load_variants",
        action="store_false",
        help="skip the catalog load, if the variants are already in",
    )
    imp.set_defaults(func=cmd_import_results, load_variants=True)

    sweep = subparsers.add_parser("sweep", help="run or resume a parameter sweep")
    sweep.add_argument("--grid", default=DEFAULT_GRID)
    sweep.add_argument("--engine", default="stub", choices=["stub", "nupack"])
    sweep.add_argument("--target-name", default="target")
    sweep.add_argument("--target-sequence", required=True)
    sweep.add_argument("--resume", type=int, default=None, metavar="RUN_ID")
    sweep.add_argument("--batch-size", type=int, default=500)
    sweep.add_argument("--note", default=None)
    sweep.add_argument("--quiet", action="store_true")
    sweep.set_defaults(func=cmd_sweep)

    rank = subparsers.add_parser("rank", help="rank by delta_g_binding and store the shortlist")
    rank.add_argument("--run-id", type=int, default=None)
    rank.add_argument("--config", default=DEFAULT_SELECTION)
    rank.add_argument("--query", default=DEFAULT_QUERY)
    rank.add_argument(
        "--strongest",
        action="store_true",
        help="rank most negative first instead, the strongest binders",
    )
    rank.add_argument("--out", default=None, help="also write the shortlist to CSV")
    rank.set_defaults(func=cmd_rank)

    trends = subparsers.add_parser("trends", help="aggregate binding over a design dimension")
    trends.add_argument(
        "--by",
        nargs="+",
        default=["truncation"],
        choices=sorted(DIMENSIONS),
        metavar="DIMENSION",
        help="one or more dimensions, see --list",
    )
    trends.add_argument("--list", action="store_true", help="show the available dimensions")
    trends.add_argument("--run-id", type=int, default=None)
    trends.add_argument("--config", default=DEFAULT_SELECTION)
    trends.add_argument(
        "--all-conditions",
        action="store_true",
        help="pool every condition instead of using the one in the config",
    )
    trends.add_argument("--min-group-size", type=int, default=1)
    trends.add_argument("--out", default=None, help="directory to write one CSV per dimension")
    trends.set_defaults(func=cmd_trends)

    explain = subparsers.add_parser("explain", help="EXPLAIN QUERY PLAN for the ranking query")
    explain.add_argument("--run-id", type=int, default=None)
    explain.add_argument("--config", default=DEFAULT_SELECTION)
    explain.add_argument("--query", default=DEFAULT_QUERY)
    explain.set_defaults(func=cmd_explain)

    runs = subparsers.add_parser("runs", help="list runs, conditions and selections")
    runs.set_defaults(func=cmd_runs)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
