-- Binding free energy aggregated over one design dimension.
--
-- Two placeholders in this file are filled by trends.py before execution:
-- {dimension_expression} and {order_by}. They are not bind parameters, because
-- a bind parameter cannot appear in GROUP BY or ORDER BY. They come from the
-- DIMENSIONS table in trends.py and never from the caller. Every value that
-- does come from the caller is bound normally.
--
-- Parameter names are listed without their leading colon, since SQLAlchemy's
-- text() scans comments for bind markers too.
--
--   run_id           run to aggregate
--   condition_id     one condition, or NULL to pool every condition in the run
--   min_group_size   drop buckets with fewer rows than this
--
-- The ddG column is the shift against the unmutated design at the same
-- truncation, so a positive mean says the changes in that bucket cost binding.

WITH scoped AS (
    SELECT
        r.variant_id,
        r.condition_id,
        r.delta_g_binding_kcal
    FROM results AS r
    WHERE r.run_id = :run_id
      AND (:condition_id IS NULL OR r.condition_id = :condition_id)
),

reference AS (
    -- One unmutated baseline per truncation per condition.
    --
    -- Driven from variants so the small set of unmutated rows is found through
    -- idx_variants_reference_lookup instead of by reading every result in the
    -- run and checking num_mutations on each one.
    SELECT
        r.condition_id,
        v.trunc_5prime,
        v.trunc_3prime,
        r.delta_g_binding_kcal AS reference_delta_g_binding_kcal
    FROM variants AS v
    JOIN results AS r
      ON r.variant_id = v.variant_id
     AND r.run_id = :run_id
     AND (:condition_id IS NULL OR r.condition_id = :condition_id)
    WHERE v.num_mutations = 0
),

annotated AS (
    SELECT
        {dimension_expression} AS bucket,
        s.delta_g_binding_kcal,
        s.delta_g_binding_kcal - ref.reference_delta_g_binding_kcal AS ddg_vs_reference_kcal
    FROM scoped AS s
    JOIN variants AS v ON v.variant_id = s.variant_id
    JOIN conditions AS c ON c.condition_id = s.condition_id
    LEFT JOIN reference AS ref
           ON ref.condition_id = s.condition_id
          AND ref.trunc_5prime = v.trunc_5prime
          AND ref.trunc_3prime = v.trunc_3prime
)

SELECT
    bucket,
    COUNT(*) AS n,
    ROUND(AVG(delta_g_binding_kcal), 3)  AS mean_delta_g_binding,
    ROUND(MIN(delta_g_binding_kcal), 3)  AS min_delta_g_binding,
    ROUND(MAX(delta_g_binding_kcal), 3)  AS max_delta_g_binding,
    -- Population standard deviation, written out because SQLite has no STDDEV.
    ROUND(
        SQRT(
            AVG(delta_g_binding_kcal * delta_g_binding_kcal)
            - AVG(delta_g_binding_kcal) * AVG(delta_g_binding_kcal)
        ),
        3
    ) AS sd_delta_g_binding,
    ROUND(AVG(ddg_vs_reference_kcal), 3) AS mean_ddg_vs_reference,
    ROUND(MAX(ddg_vs_reference_kcal), 3) AS max_ddg_vs_reference
FROM annotated
WHERE bucket IS NOT NULL
GROUP BY bucket
HAVING COUNT(*) >= :min_group_size
ORDER BY {order_by};
