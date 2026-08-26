-- Rank designs by binding free energy, most negative first.
--
-- This is the mirror of sql/select_candidates.sql. Sort direction cannot be a
-- bind parameter, so the two directions are two files rather than one file with
-- a flag. Everything else, including the parameter list, is identical, and a
-- test checks the two stay in step.
--
-- Use this one to find the strongest binders. Use select_candidates.sql to find
-- the designs whose binding is most disrupted, which is the usual question on a
-- mutation scan.
--
-- Parameter names are listed below without their leading colon on purpose.
-- SQLAlchemy's text() scans the whole file for bind markers and does not skip
-- comments, so writing a marker up here would register an extra parameter and
-- the execute would fail on a binding count mismatch. Keep it that way.
--
--   run_id                   run to select from
--   condition_id             condition to rank at, resolved in Python
--   min_delta_g_binding      lower bound on binding free energy
--   max_delta_g_binding      upper bound on binding free energy
--   min_ddg                  lower bound on the shift against the reference
--   gc_min, gc_max           GC content window
--   trunc_total_min          fewest bases trimmed
--   trunc_total_max          most bases trimmed
--   num_mutations_min        fewest mutations
--   num_mutations_max        most mutations
--   per_class_cap            most rows to take from any one variant class
--   n_final                  size of the shortlist
--
-- Shape: pull the numbers at one condition, attach the matching unmutated
-- reference so each row carries a ddG, filter, then rank within class.

WITH scored AS (
    -- Hits idx_results_ranking, which also covers the ORDER BY at the bottom.
    SELECT
        r.variant_id,
        r.delta_g_binding_kcal,
        r.cdna_dg_kcal,
        r.complex_dg_kcal
    FROM results AS r
    WHERE r.run_id = :run_id
      AND r.condition_id = :condition_id
),

reference AS (
    -- The unmutated design at each truncation, which is the right baseline for
    -- a mutation effect. Comparing everything to the untruncated wild type
    -- would fold the cost of truncation into every ddG and swamp the mutation
    -- signal, since a scan usually trims several bases.
    --
    -- Driven from variants rather than from the scored CTE, which states the
    -- intent directly: find the handful of unmutated designs, then fetch their
    -- results. SQLite reorders both forms to the same plan on this data, so the
    -- choice is about readability rather than speed.
    SELECT
        v.trunc_5prime,
        v.trunc_3prime,
        r.delta_g_binding_kcal AS reference_delta_g_binding_kcal
    FROM variants AS v
    JOIN results AS r
      ON r.variant_id = v.variant_id
     AND r.run_id = :run_id
     AND r.condition_id = :condition_id
    WHERE v.num_mutations = 0
),

annotated AS (
    SELECT
        v.variant_id,
        v.name,
        v.sequence,
        v.variant_class,
        v.mutation_signature,
        v.length_nt,
        v.gc_content,
        v.trunc_5prime,
        v.trunc_3prime,
        v.trunc_total,
        v.num_mutations,
        v.mutation_type_1,
        v.mutation_type_2,
        v.num_transitions,
        v.num_transversions,
        s.cdna_dg_kcal,
        s.complex_dg_kcal,
        s.delta_g_binding_kcal,
        ref.reference_delta_g_binding_kcal,
        -- Positive ddG means the change cost binding relative to the unmutated
        -- design at the same truncation. LEFT JOIN so a scan with no reference
        -- row still ranks, it just has no ddG.
        s.delta_g_binding_kcal - ref.reference_delta_g_binding_kcal AS ddg_vs_reference_kcal
    FROM scored AS s
    JOIN variants AS v ON v.variant_id = s.variant_id
    LEFT JOIN reference AS ref
           ON ref.trunc_5prime = v.trunc_5prime
          AND ref.trunc_3prime = v.trunc_3prime
),

filtered AS (
    SELECT *
    FROM annotated
    WHERE delta_g_binding_kcal BETWEEN :min_delta_g_binding AND :max_delta_g_binding
      AND gc_content BETWEEN :gc_min AND :gc_max
      AND trunc_total BETWEEN :trunc_total_min AND :trunc_total_max
      AND num_mutations BETWEEN :num_mutations_min AND :num_mutations_max
      AND (ddg_vs_reference_kcal IS NULL OR ddg_vs_reference_kcal >= :min_ddg)
),

ranked AS (
    SELECT
        filtered.*,
        ROW_NUMBER() OVER (
            PARTITION BY variant_class
            ORDER BY delta_g_binding_kcal ASC, variant_id ASC
        ) AS rank_in_class
    FROM filtered
)

SELECT
    variant_id,
    name,
    sequence,
    variant_class,
    mutation_signature,
    length_nt,
    trunc_5prime,
    trunc_3prime,
    trunc_total,
    num_mutations,
    num_transitions,
    num_transversions,
    ROUND(gc_content, 4)                  AS gc_content,
    ROUND(cdna_dg_kcal, 3)                AS cdna_dg_kcal,
    ROUND(complex_dg_kcal, 3)             AS complex_dg_kcal,
    ROUND(delta_g_binding_kcal, 3)        AS delta_g_binding_kcal,
    ROUND(reference_delta_g_binding_kcal, 3) AS reference_delta_g_binding_kcal,
    ROUND(ddg_vs_reference_kcal, 3)       AS ddg_vs_reference_kcal,
    rank_in_class
FROM ranked
WHERE rank_in_class <= :per_class_cap
ORDER BY delta_g_binding_kcal ASC, variant_id ASC
LIMIT :n_final;
