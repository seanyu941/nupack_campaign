"""Run the selection query and record what came out of it.

The query itself lives in ``sql/select_candidates.sql``. This module is the thin
layer around it: read the config, turn the named condition into a condition_id,
bind the parameters, and write the shortlist back to the database so the result
can be looked at again later without rerunning anything.

Parameters are always bound, never formatted into the string. Apart from the
injection question, bound parameters let SQLite reuse the prepared statement
when the same query is run across several thresholds.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
import yaml
from sqlalchemy import func, select, text
from sqlalchemy.orm import Session

from .models import Condition, Result, Selection, SelectionMember

DEFAULT_QUERY_PATH = Path("sql/select_candidates.sql")
STRONGEST_QUERY_PATH = Path("sql/select_strongest_binders.sql")
DEFAULT_CONFIG_PATH = Path("config/selection.yaml")

# Condition values come from YAML and land in SQLite as REAL, so comparing them
# back is an exact float comparison. Resolving through a tolerance avoids the
# usual 0.1 + 0.2 surprise.
FLOAT_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class SelectionCriteria:
    """Flat view of config/selection.yaml, one field per bind parameter."""

    temperature_c: float
    na_molar: float
    mg_molar: float
    material: str
    ensemble: str
    min_delta_g_binding: float
    max_delta_g_binding: float
    min_ddg: float
    gc_min: float
    gc_max: float
    trunc_total_min: int
    trunc_total_max: int
    num_mutations_min: int
    num_mutations_max: int
    per_class_cap: int
    n_final: int

    @classmethod
    def from_yaml(cls, path: str | Path = DEFAULT_CONFIG_PATH) -> SelectionCriteria:
        raw = yaml.safe_load(Path(path).read_text())
        condition = raw["condition"]
        thresholds = raw["thresholds"]
        limits = raw["limits"]

        criteria = cls(
            temperature_c=float(condition["temperature_c"]),
            na_molar=float(condition["na_molar"]),
            mg_molar=float(condition["mg_molar"]),
            material=str(condition["material"]),
            ensemble=str(condition["ensemble"]),
            min_delta_g_binding=float(thresholds["min_delta_g_binding"]),
            max_delta_g_binding=float(thresholds["max_delta_g_binding"]),
            min_ddg=float(thresholds["min_ddg"]),
            gc_min=float(thresholds["gc_min"]),
            gc_max=float(thresholds["gc_max"]),
            trunc_total_min=int(thresholds["trunc_total_min"]),
            trunc_total_max=int(thresholds["trunc_total_max"]),
            num_mutations_min=int(thresholds["num_mutations_min"]),
            num_mutations_max=int(thresholds["num_mutations_max"]),
            per_class_cap=int(limits["per_class_cap"]),
            n_final=int(limits["n_final"]),
        )
        criteria.validate()
        return criteria

    def validate(self) -> None:
        if self.min_delta_g_binding > self.max_delta_g_binding:
            raise ValueError("min_delta_g_binding is above max_delta_g_binding, the band is empty")
        if self.gc_min > self.gc_max:
            raise ValueError("gc_min is above gc_max")
        if self.trunc_total_min > self.trunc_total_max:
            raise ValueError("truncation bounds are the wrong way round")
        if self.num_mutations_min > self.num_mutations_max:
            raise ValueError("mutation count bounds are the wrong way round")
        if self.per_class_cap < 1 or self.n_final < 1:
            raise ValueError("per_class_cap and n_final should both be at least 1")


class ConditionNotFound(LookupError):
    pass


def resolve_condition_id(session: Session, criteria: SelectionCriteria) -> int:
    """Find the condition row named in the config.

    Raising here rather than returning an empty shortlist is deliberate. An
    empty result from a condition that was never simulated looks exactly like a
    condition where nothing passed the filters, and those are very different
    problems.
    """
    condition_id = session.execute(
        select(Condition.condition_id).where(
            func.abs(Condition.temperature_c - criteria.temperature_c) < FLOAT_TOLERANCE,
            func.abs(Condition.na_molar - criteria.na_molar) < FLOAT_TOLERANCE,
            func.abs(Condition.mg_molar - criteria.mg_molar) < FLOAT_TOLERANCE,
            Condition.material == criteria.material,
            Condition.ensemble == criteria.ensemble,
        )
    ).scalar_one_or_none()

    if condition_id is None:
        available = session.execute(
            select(Condition.temperature_c, Condition.na_molar, Condition.mg_molar)
        ).all()
        raise ConditionNotFound(
            f"No condition at {criteria.temperature_c} C, Na={criteria.na_molar} M, "
            f"Mg={criteria.mg_molar} M, {criteria.material}/{criteria.ensemble}. "
            f"The database has: {sorted(set(available))[:6]}"
        )
    return int(condition_id)


def load_query(path: str | Path = DEFAULT_QUERY_PATH) -> str:
    return Path(path).read_text()


def bind_parameters(criteria: SelectionCriteria, run_id: int, condition_id: int) -> dict:
    return {
        "run_id": run_id,
        "condition_id": condition_id,
        "min_delta_g_binding": criteria.min_delta_g_binding,
        "max_delta_g_binding": criteria.max_delta_g_binding,
        "min_ddg": criteria.min_ddg,
        "gc_min": criteria.gc_min,
        "gc_max": criteria.gc_max,
        "trunc_total_min": criteria.trunc_total_min,
        "trunc_total_max": criteria.trunc_total_max,
        "num_mutations_min": criteria.num_mutations_min,
        "num_mutations_max": criteria.num_mutations_max,
        "per_class_cap": criteria.per_class_cap,
        "n_final": criteria.n_final,
    }


def select_candidates(
    session: Session,
    run_id: int,
    criteria: SelectionCriteria,
    query_path: str | Path = DEFAULT_QUERY_PATH,
) -> pd.DataFrame:
    """Execute the shortlist query and return it as a DataFrame.

    Defaults to ``sql/select_candidates.sql``, which ranks most positive
    delta_g_binding first. Pass ``STRONGEST_QUERY_PATH`` for the opposite end.
    """
    condition_id = resolve_condition_id(session, criteria)
    sql = load_query(query_path)
    params = bind_parameters(criteria, run_id, condition_id)

    result = session.execute(text(sql), params)
    # Column names are read before the rows are drained, so a shortlist that
    # comes back empty still has the right columns. Otherwise callers have to
    # check for an empty frame before touching any column by name.
    columns = list(result.keys())
    return pd.DataFrame(result.mappings().all(), columns=columns)


def count_candidates(session: Session, run_id: int, condition_id: int) -> int:
    """How many variants were in play at that condition before filtering."""
    return int(
        session.execute(
            select(func.count(func.distinct(Result.variant_id))).where(
                Result.run_id == run_id, Result.condition_id == condition_id
            )
        ).scalar_one()
    )


def persist_selection(
    session: Session,
    run_id: int,
    criteria: SelectionCriteria,
    shortlist: pd.DataFrame,
    query_path: str | Path = DEFAULT_QUERY_PATH,
) -> Selection:
    """Write the shortlist and the criteria that produced it.

    The criteria are stored as a snapshot rather than a path to the config file,
    because the config file will change and the shortlist should still explain
    itself six months later.
    """
    condition_id = resolve_condition_id(session, criteria)
    sql = load_query(query_path)

    selection = Selection(
        run_id=run_id,
        criteria_json=json.dumps(asdict(criteria), sort_keys=True),
        sql_sha256=hashlib.sha256(sql.encode()).hexdigest(),
        n_candidates_in=count_candidates(session, run_id, condition_id),
        n_selected=len(shortlist),
    )
    session.add(selection)
    session.flush()

    for position, row in enumerate(shortlist.to_dict("records"), start=1):
        ddg = row.get("ddg_vs_reference_kcal")
        session.add(
            SelectionMember(
                selection_id=selection.selection_id,
                variant_id=int(row["variant_id"]),
                rank=position,
                delta_g_binding_kcal=float(row["delta_g_binding_kcal"]),
                ddg_vs_reference_kcal=None if ddg is None or pd.isna(ddg) else float(ddg),
            )
        )

    session.commit()
    return selection


def explain_query_plan(
    session: Session,
    run_id: int,
    criteria: SelectionCriteria,
    query_path: str | Path = DEFAULT_QUERY_PATH,
) -> str:
    """EXPLAIN QUERY PLAN output, used to check the indexes are doing something.

    SQLite only. Worth rerunning after any change to the query or the indexes.
    """
    condition_id = resolve_condition_id(session, criteria)
    sql = load_query(query_path).rstrip().rstrip(";")
    params = bind_parameters(criteria, run_id, condition_id)

    rows = session.execute(text(f"EXPLAIN QUERY PLAN {sql}"), params).all()
    return "\n".join(str(row[3]) for row in rows)
