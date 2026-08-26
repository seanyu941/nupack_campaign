"""Trend queries: how binding moves with one design dimension at a time.

The aggregate is always the same, only the grouping expression changes, so the
query is assembled from ``sql/trends/trend_by_dimension.sql`` with the GROUP BY
expression substituted in.

That substitution is string formatting into SQL, which is normally the thing to
avoid. It is safe here only because the expression never comes from the caller:
the dimension name is looked up in ``DIMENSIONS`` below and the SQL comes from
this file. Anything not in the table raises. Every value the caller does supply
is still a bind parameter. If you add a dimension, write the expression here
rather than accepting one from a config file or an argument.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sqlalchemy import text
from sqlalchemy.orm import Session

DEFAULT_TREND_QUERY = Path("sql/trends/trend_by_dimension.sql")


@dataclass(frozen=True, slots=True)
class Dimension:
    """One grouping axis.

    ``expression`` is SQL evaluated against the joined variants and results
    rows. ``label`` is what the output column is called. ``order_by`` is how the
    buckets are sorted, which for a binned dimension is the bucket value rather
    than the aggregate.
    """

    expression: str
    label: str
    description: str
    order_by: str = "bucket"


DIMENSIONS: dict[str, Dimension] = {
    "truncation": Dimension(
        expression="v.trunc_total",
        label="bases trimmed",
        description="total bases removed from either end",
    ),
    "truncation_5prime": Dimension(
        expression="v.trunc_5prime",
        label="5' bases trimmed",
        description="bases removed from the 5' end",
    ),
    "truncation_3prime": Dimension(
        expression="v.trunc_3prime",
        label="3' bases trimmed",
        description="bases removed from the 3' end",
    ),
    "truncation_pair": Dimension(
        expression="v.trunc_5prime || '/' || v.trunc_3prime",
        label="5'/3' trimmed",
        description="the two truncation ends as a pair",
    ),
    "gc": Dimension(
        # 0.02 wide buckets. Raw gc_content has one distinct value per base
        # count per length, which is too many rows to read as a trend.
        expression="ROUND(v.gc_content / 0.02) * 0.02",
        label="GC content bin",
        description="GC content in bins of 0.02",
    ),
    "length": Dimension(
        expression="v.length_nt",
        label="length (nt)",
        description="design length after truncation",
    ),
    "num_mutations": Dimension(
        expression="v.num_mutations",
        label="mutations",
        description="how many positions were substituted",
    ),
    "mutation_type": Dimension(
        # Ordered so that transition/transversion and transversion/transition
        # land in the same bucket, since the slots are not meaningfully ordered.
        expression=(
            "CASE WHEN v.mutation_type_1 <= v.mutation_type_2 "
            "THEN v.mutation_type_1 || '+' || v.mutation_type_2 "
            "ELSE v.mutation_type_2 || '+' || v.mutation_type_1 END"
        ),
        label="mutation types",
        description="the pair of substitution types, order independent",
    ),
    "transitions": Dimension(
        expression="v.num_transitions",
        label="transitions",
        description="how many substitutions were transitions",
    ),
    "transversions": Dimension(
        expression="v.num_transversions",
        label="transversions",
        description="how many substitutions were transversions",
    ),
    "position": Dimension(
        expression="v.position_1",
        label="first mutated position",
        description="position of the first substitution",
    ),
    "position_2": Dimension(
        expression="v.position_2",
        label="second mutated position",
        description="position of the second substitution",
    ),
    "variant_class": Dimension(
        expression="v.variant_class",
        label="variant class",
        description="mutation and truncation bucket",
    ),
    "original_base": Dimension(
        expression="v.original_base_1",
        label="original base",
        description="base that was replaced at the first position",
    ),
    "substitution": Dimension(
        expression="v.original_base_1 || '>' || v.mutated_base_1",
        label="substitution",
        description="the first substitution, such as G>A",
    ),
    "paired_bases": Dimension(
        expression="a.n_complementary",
        label="paired bases",
        description="bases that actually pair with the target, truncation and mutation combined",
    ),
    "mismatches": Dimension(
        expression="a.n_mismatch",
        label="mismatches",
        description="bases in the binding window that do not pair",
    ),
    "binding_offset": Dimension(
        expression="a.target_offset",
        label="target offset",
        description="where the binding window starts on the target",
    ),
    "match_fraction": Dimension(
        expression="ROUND(a.match_fraction / 0.02) * 0.02",
        label="paired fraction bin",
        description="fraction of the design that pairs, in bins of 0.02",
    ),
    "temperature": Dimension(
        expression="c.temperature_c",
        label="temperature (C)",
        description="sweep temperature, for trends across conditions",
    ),
}


class UnknownDimension(KeyError):
    pass


def resolve(dimension: str) -> Dimension:
    try:
        return DIMENSIONS[dimension]
    except KeyError:
        options = ", ".join(sorted(DIMENSIONS))
        raise UnknownDimension(
            f"{dimension!r} is not a known dimension. Pick one of: {options}"
        ) from None


def describe_dimensions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"dimension": name, "groups by": dim.description}
            for name, dim in sorted(DIMENSIONS.items())
        ]
    )


def load_trend_query(path: str | Path = DEFAULT_TREND_QUERY) -> str:
    return Path(path).read_text()


def trend_by(
    session: Session,
    run_id: int,
    dimension: str,
    condition_id: int | None = None,
    min_group_size: int = 1,
    query_path: str | Path = DEFAULT_TREND_QUERY,
) -> pd.DataFrame:
    """Aggregate binding free energy over one dimension.

    Pass ``condition_id`` to look at a single condition. Leave it as None to
    pool every condition in the run, which is what the temperature dimension
    needs.
    """
    dim = resolve(dimension)
    sql = load_trend_query(query_path).format(
        dimension_expression=dim.expression,
        order_by=dim.order_by,
    )

    params = {
        "run_id": run_id,
        "condition_id": condition_id,
        "min_group_size": min_group_size,
    }

    result = session.execute(text(sql), params)
    columns = list(result.keys())
    frame = pd.DataFrame(result.mappings().all(), columns=columns)
    return frame.rename(columns={"bucket": dim.label})
