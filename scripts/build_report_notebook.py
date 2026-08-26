"""Write notebooks/01_campaign_report.ipynb.

The notebook is generated rather than hand edited so it stays free of execution
counts and stored outputs, which are the two things that make notebook diffs
unreadable. Regenerate after changing the cells here:

    python scripts/build_report_notebook.py
"""

from __future__ import annotations

import json
from pathlib import Path

OUT_PATH = Path("notebooks/01_campaign_report.ipynb")

MARKDOWN_INTRO = """\
# Campaign report

Reads the scan out of the database and looks at it. Nothing here computes
anything or writes anything, so it can be rerun against any database without
touching the results.

The order the stages run in is:

1. `make import` loads the scan CSV and unpivots its energy columns, or
   `make sweep` simulates them instead
2. `make rank` runs `sql/select_candidates.sql` and stores the shortlist
3. `make trends` aggregates over the design dimensions
4. this notebook reads all of it back

If `data/campaign.db` is not there yet, run `make demo` from the repo root
first.\
"""

SETUP_CODE = """\
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sqlalchemy import select, text

REPO_ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(REPO_ROOT / "src"))

from nupack_campaign.db import make_engine, session_scope
from nupack_campaign.models import Selection
from nupack_campaign.selection import (
    SelectionCriteria,
    resolve_condition_id,
    select_candidates,
)
from nupack_campaign.trends import DIMENSIONS, trend_by

engine = make_engine(f"sqlite:///{REPO_ROOT / 'data' / 'campaign.db'}")
SQL = REPO_ROOT / "sql"
pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 30)\
"""

MARKDOWN_RUNS = """\
## What is in the database

One row per sweep. `grid_hash` covers the contents of `config/sweep_grid.yaml`
and `git_sha` the commit, so a result set can be traced back to the code and
the grid that produced it.\
"""

RUNS_CODE = """\
runs = pd.read_sql(
    text(
        \"\"\"
        SELECT r.run_id, r.status, r.engine_name, r.engine_version,
               SUBSTR(r.git_sha, 1, 8) AS git_sha,
               SUBSTR(r.grid_hash, 1, 8) AS grid_hash,
               r.started_at, r.finished_at,
               COUNT(res.result_id) AS n_results,
               ROUND(SUM(res.compute_seconds), 1) AS compute_seconds
        FROM runs r
        LEFT JOIN results res ON res.run_id = r.run_id
        GROUP BY r.run_id
        ORDER BY r.run_id
        \"\"\"
    ),
    engine,
)
runs\
"""

MARKDOWN_COVERAGE = """\
## Scan coverage

The scan is a grid of its own: truncation from either end crossed with every
double substitution. This pivot is a sanity check that no corner of it is
missing.\
"""

COVERAGE_CODE = """\
RUN_ID = int(runs.loc[runs.status == "complete", "run_id"].max())

coverage = pd.read_sql(
    text(
        \"\"\"
        SELECT v.trunc_5prime, v.trunc_3prime,
               COUNT(*) AS n_variants,
               SUM(CASE WHEN v.num_mutations = 0 THEN 1 ELSE 0 END) AS n_reference
        FROM variants v
        GROUP BY v.trunc_5prime, v.trunc_3prime
        ORDER BY v.trunc_5prime, v.trunc_3prime
        \"\"\"
    ),
    engine,
)

print(f"{coverage.n_variants.sum()} variants, {coverage.n_reference.sum()} unmutated references")
coverage.pivot_table(index="trunc_5prime", columns="trunc_3prime", values="n_variants")\
"""

MARKDOWN_RANKING = """\
## The ranking

`sql/select_candidates.sql` sorts by `delta_g_binding` descending, most positive
first, so the top of the list is where binding is most disrupted. The mirror
query sorts the other way for the strongest binders.

`ddg_vs_reference_kcal` is measured against the unmutated design at the same
truncation, not against the untruncated wild type. Truncation costs about
1.1 kcal/mol per base on this scan, which would otherwise swamp the mutation
signal in every ddG.\
"""

RANKING_CODE = """\
criteria = SelectionCriteria.from_yaml(REPO_ROOT / "config" / "selection.yaml")

with session_scope(engine) as session:
    condition_id = resolve_condition_id(session, criteria)
    disrupted = select_candidates(
        session, RUN_ID, criteria, query_path=SQL / "select_candidates.sql"
    )
    strongest = select_candidates(
        session, RUN_ID, criteria, query_path=SQL / "select_strongest_binders.sql"
    )

columns = [
    "name",
    "variant_class",
    "mutation_signature",
    "trunc_5prime",
    "trunc_3prime",
    "gc_content",
    "delta_g_binding_kcal",
    "ddg_vs_reference_kcal",
]

print("Most disrupted binding:")
display(disrupted[columns].head(10))
print("\\nStrongest binding:")
display(strongest[columns].head(10))\
"""

MARKDOWN_TRENDS = """\
## Trends

One aggregate per design dimension, all going through the same template in
`sql/trends/trend_by_dimension.sql`. `DIMENSIONS` in `trends.py` lists what can
be grouped on.\
"""

TRENDS_CODE = """\
with session_scope(engine) as session:
    trends = {
        name: trend_by(
            session,
            RUN_ID,
            name,
            condition_id=condition_id,
            query_path=SQL / "trends" / "trend_by_dimension.sql",
        )
        for name in ("truncation", "paired_bases", "gc", "position", "substitution")
    }

for name, frame in trends.items():
    print(f"\\n{name}: {DIMENSIONS[name].description}")
    display(frame.head(12))\
"""

MARKDOWN_TREND_PLOT = """\
## Trends, drawn

Left: binding tracks the number of bases that actually pair, which is where
truncation and mutation both end up. Middle: a mutation near either end costs
less than one in the middle, the expected shape for a terminal mismatch.\
"""

TREND_PLOT_CODE = """\
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

trunc = trends["paired_bases"]
axes[0].errorbar(
    trunc.iloc[:, 0], trunc.mean_delta_g_binding, yerr=trunc.sd_delta_g_binding,
    marker="o", capsize=3,
)
axes[0].set_xlabel("bases paired with the target")
axes[0].set_ylabel("mean delta_g_binding (kcal/mol)")
axes[0].set_title("Pairing drives binding")

position = trends["position"]
axes[1].plot(position.iloc[:, 0], position.mean_ddg_vs_reference, marker="o", color="tab:orange")
axes[1].set_xlabel("first mutated position")
axes[1].set_ylabel("mean ddG vs reference (kcal/mol)")
axes[1].set_title("Where a mutation hurts most")

substitution = trends["substitution"].sort_values("mean_ddg_vs_reference", ascending=False)
axes[2].bar(substitution.iloc[:, 0], substitution.mean_ddg_vs_reference, color="tab:green")
axes[2].set_xlabel("substitution")
axes[2].set_ylabel("mean ddG vs reference (kcal/mol)")
axes[2].set_title("Which substitution hurts most")
axes[2].tick_params(axis="x", rotation=45)

for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()\
"""

MARKDOWN_LANDSCAPE = """\
## The whole landscape

Every variant at the ranking condition, coloured by how far it was truncated.
The shortlist sits at the right hand edge, which is what ranking most positive
first picks out.\
"""

LANDSCAPE_CODE = """\
landscape = pd.read_sql(
    text(
        \"\"\"
        SELECT v.trunc_total, v.gc_content, r.delta_g_binding_kcal
        FROM results r
        JOIN variants v ON v.variant_id = r.variant_id
        WHERE r.run_id = :run_id AND r.condition_id = :condition_id
        \"\"\"
    ),
    engine,
    params={"run_id": RUN_ID, "condition_id": condition_id},
)

fig, ax = plt.subplots(figsize=(10, 5))
points = ax.scatter(
    landscape.delta_g_binding_kcal,
    landscape.gc_content,
    c=landscape.trunc_total,
    s=6,
    alpha=0.4,
    cmap="viridis",
)
ax.scatter(
    disrupted.delta_g_binding_kcal,
    disrupted.gc_content,
    s=90,
    facecolor="none",
    edgecolor="crimson",
    linewidth=1.4,
    label="shortlist",
    zorder=3,
)
fig.colorbar(points, ax=ax, label="bases trimmed")
ax.set_xlabel("delta_g_binding (kcal/mol)")
ax.set_ylabel("GC content")
ax.set_title(f"Run {RUN_ID}: {len(landscape)} variants, shortlist ringed")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()\
"""

MARKDOWN_PROVENANCE = """\
## Stored selections

Every run of `make rank` appends a row here with the criteria it used, so an old
shortlist can be explained without rerunning it or digging through git history
for what the config said at the time.\
"""

PROVENANCE_CODE = """\
import json

with session_scope(engine) as session:
    records = session.execute(select(Selection).order_by(Selection.selection_id)).scalars().all()
    rows = [
        {
            "selection_id": record.selection_id,
            "run_id": record.run_id,
            "in": record.n_candidates_in,
            "out": record.n_selected,
            "sql_sha256": record.sql_sha256[:8],
            **{
                key: value
                for key, value in json.loads(record.criteria_json).items()
                if key.startswith(("min_", "max_", "gc_", "trunc_", "num_", "per_", "n_"))
            },
        }
        for record in records
    ]

pd.DataFrame(rows)\
"""


def markdown_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def build_notebook() -> dict:
    cells = [
        markdown_cell(MARKDOWN_INTRO),
        code_cell(SETUP_CODE),
        markdown_cell(MARKDOWN_RUNS),
        code_cell(RUNS_CODE),
        markdown_cell(MARKDOWN_COVERAGE),
        code_cell(COVERAGE_CODE),
        markdown_cell(MARKDOWN_RANKING),
        code_cell(RANKING_CODE),
        markdown_cell(MARKDOWN_TRENDS),
        code_cell(TRENDS_CODE),
        markdown_cell(MARKDOWN_TREND_PLOT),
        code_cell(TREND_PLOT_CODE),
        markdown_cell(MARKDOWN_LANDSCAPE),
        code_cell(LANDSCAPE_CODE),
        markdown_cell(MARKDOWN_PROVENANCE),
        code_cell(PROVENANCE_CODE),
    ]

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.10"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(build_notebook(), indent=1) + "\n")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
