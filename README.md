# NUPACK Campaign

Takes a DNA design scan, stores every free energy result in SQLite, and answers
two kinds of question about it with SQL:

- **Ranking.** Which designs sit where on binding free energy, sorted most
  positive first so the top of the list is where binding is most disrupted.
- **Trends.** How binding moves with one design dimension at a time: truncation,
  GC content, mutation position, substitution type, temperature.

Built around a 25 nt cDNA double mutation and truncation scan: 54,660 designs
across 2 temperatures, 109,320 result rows.

```
[load-variants] 54660 inserted (double_mutation 2700, double_mutation_truncated 51930,
                                reference 1, truncation_only 29)
[import] run 1: 109320 rows across 2 conditions (25C, 37C)
[import] recovered target dG at 25C: -27.4946 kcal/mol
[import] recovered target dG at 37C: -21.3736 kcal/mol
[rank] selection 1: 54660 candidates -> 16 selected, ranked by delta_g_binding most positive first
```

## Quick start

```bash
pip install -e ".[dev]"
make demo
```

`make demo` imports the scan CSV, ranks it, and prints the trends. It takes
about 30 seconds and needs no NUPACK licence, because the CSV already has the
free energies in it.

## Two ways in

**If your scan CSV already has free energies**, import it. Nothing is
recomputed:

```bash
python -m nupack_campaign.cli import-results \
    --catalog data/scan_25nt.csv \
    --target-sequence TTGCTGCGCGCTGAAGCGCGGAAGC
```

**If it does not**, load the designs and simulate:

```bash
python -m nupack_campaign.cli load-variants --catalog designs.csv
python -m nupack_campaign.cli sweep --engine nupack --target-sequence ...
```

Either way the result lands in the same tables, so an imported run and a
simulated run sit side by side and every query works on both.

## Ranking

```bash
make rank                    # most positive delta_g_binding first
make strongest               # the mirror, most negative first
```

```
          name             variant_class mutation_signature  trunc_5  trunc_3  gc     delta_g_binding  ddg_vs_reference
  t5_4_T5G,G9C double_mutation_truncated            T5G,G9C        5        4  0.8125          -1.512             9.354
  t5_4_T5G,G9T double_mutation_truncated            T5G,G9T        5        4  0.7500          -1.791             9.075
 t5_4_C7A,C10A double_mutation_truncated           C7A,C10A        5        4  0.6250          -2.231             8.635
 t0_0_T3G,T10C           double_mutation           T3G,T10C        0        0  0.7600         -10.218            10.961
       t5_4_wt           truncation_only                           5        4  0.7500         -10.866             0.000
       t0_0_wt                 reference                           0        0  0.6800         -21.179             0.000
```

Sort direction is `DESC` in `sql/select_candidates.sql` and that is deliberate.
Binding free energy is negative for a design that binds, so most positive means
most disrupted. There is a test named
`test_ranking_is_most_positive_first` guarding it, because a passing glance at
`ORDER BY ... DESC` looks like something to tidy up.

Direction cannot be a bind parameter, so the opposite ranking is a second file,
`sql/select_strongest_binders.sql`. A test checks the two take an identical
parameter list.

### ddG is measured against the matching truncation

Each truncation combination has its own unmutated reference row, 30 of them on
this scan. `ddg_vs_reference_kcal` compares a mutant to the reference at the
*same* truncation rather than to the untruncated wild type. Truncation costs
about 1.1 kcal/mol per base here, so using one global baseline would fold that
into every ddG and bury the mutation signal underneath it.

## Trends

```bash
make trends                                  # the default set
python -m nupack_campaign.cli trends --by truncation gc position
python -m nupack_campaign.cli trends --list  # what can be grouped on
```

```
 bases trimmed    n  mean_delta_g_binding  sd  mean_ddg_vs_reference  max_ddg_vs_reference
             0 2701               -16.464 1.946                  4.715                10.961
             1 4970               -15.707 1.877                  5.090                11.375
             2 6834               -14.860 1.801                  5.350                11.325
             3 8320               -13.745 1.724                  5.332                11.070
             ...
             9 1081                -6.876 1.873                  3.990                 9.354
```

Seventeen dimensions are available: `truncation`, `truncation_5prime`,
`truncation_3prime`, `truncation_pair`, `gc`, `length`, `num_mutations`,
`mutation_type`, `transitions`, `transversions`, `position`, `position_2`,
`substitution`, `original_base`, `variant_class`, `temperature`.

All of them go through one template, `sql/trends/trend_by_dimension.sql`, with
the GROUP BY expression substituted in. That substitution is string formatting
into SQL, which is normally the thing to avoid. It is safe here only because the
expression never comes from the caller: the dimension name is looked up in the
`DIMENSIONS` table in `trends.py` and anything not in the table raises. Every
value the caller does supply is still a bind parameter. A test checks each
declared expression touches only the two table aliases the template defines.

A few things the scan shows, as a check that the plumbing produces real signal:

- Truncation costs roughly 1.1 kcal/mol per base, near linear across 0 to 9.
- Mutations bite hardest at moderate truncation. Mean ddG peaks at 2 to 3 bases
  trimmed, then falls once there is little duplex left to disrupt.
- A mutation at position 0 costs 3.5 kcal/mol on average, rising to 6.0 by
  position 10. Terminal mismatches are cheap, internal ones are not.
- Transversion pairs disrupt slightly more than transition pairs, 5.14 against
  4.83 kcal/mol mean ddG.

## Schema

Eight tables. The DDL in `sql/schema.sql` and `sql/indexes.sql` is generated
from `src/nupack_campaign/models.py` by `scripts/dump_schema_sql.py`, and a test
fails if the two drift.

| table | what it holds |
| --- | --- |
| `variants` | the design catalog with its scan metadata, written once |
| `targets` | the strand the designs are meant to bind |
| `conditions` | one row per temperature and buffer combination |
| `runs` | one row per sweep or import, with git sha and grid hash |
| `target_energies` | free energy of a target on its own, per run and condition |
| `results` | one row per (run, variant, target, condition) |
| `selections` | one row per execution of a ranking query |
| `selection_members` | the shortlisted variants and what they were picked on |

### Why the wide CSV gets unpivoted

A scan CSV stores results wide, one group of columns per temperature:

```
cdna_dg_25C, complex_dg_25C, delta_g_binding_25C, cdna_dg_37C, complex_dg_37C, delta_g_binding_37C
```

Adding a third temperature means adding three more columns and touching every
query that mentions them. Asking how binding moves with temperature means
unpivoting by hand each time. `importer.py` does that once on the way in, so a
temperature becomes rows rather than columns and the trend query gets it for
free.

### Why the target energy has its own table

Binding free energy is a three-body quantity:

```
delta_g_binding = complex_dg - cdna_dg - target_dg
```

The target term does not depend on the variant, so it is one row per (run,
target, condition) rather than a copy on every result. On this scan that is 2
values instead of 109,320 duplicates of the same two numbers.

The CSV does not carry that term, so the importer backs it out from the other
three and checks it is actually constant across the file, within 0.01 kcal/mol.
On the 25 nt scan it recovers -27.4946 at 25C and -21.3736 at 37C with a spread
of 0.0005, which is the file's own rounding. If the three columns ever disagree
with each other, the spread blows up and the import stops.

### Indexes

`make explain`:

```
MATERIALIZE reference
SEARCH r USING INDEX idx_results_ranking (run_id=? AND condition_id=?)
SEARCH r USING INDEX idx_results_ranking (run_id=? AND condition_id=? AND delta_g_binding_kcal>? AND delta_g_binding_kcal<?)
BLOOM FILTER ON ref (trunc_5prime=? AND trunc_3prime=?)
SEARCH ref USING AUTOMATIC COVERING INDEX (trunc_5prime=? AND trunc_3prime=?) LEFT-JOIN
```

`idx_results_ranking` is on `(run_id, condition_id, delta_g_binding_kcal)`, so
it covers the filter and the sort. No scan of `results`, which at 109,320 rows
is the only table large enough to care about.

Two things had to change to get that plan:

1. The `reference` CTE originally drove from `results` and checked
   `num_mutations` on each row, reading all 54,660 to find 30. It now drives
   from `variants` through `idx_variants_reference_lookup` and fetches only
   those 30 results.
2. `ANALYZE` runs after every bulk load. Without table statistics SQLite picks
   the wrong join order for the reference join and falls back to scanning the
   materialised CTE. With them it builds an automatic covering index and a bloom
   filter. The trend query goes from 0.21s to 0.09s.

Current timings on the full scan: ranking 0.36s, trends 0.09s, database 26 MB.

## Engines

The sweep runner only knows about the `ThermoEngine` protocol, so there are two
implementations behind it.

`NupackEngine` is the real one, written against the NUPACK 4 API. NUPACK is not
on PyPI and needs a licence, so it is not a dependency of this package. Install
it separately from nupack.org and run `make sweep ENGINE=nupack`.

`StubEngine` is a deterministic stand-in. Values come from a hash of the
sequence and condition, so the same input always gives the same output. It
models the duplex as both strands plus a hybridisation term that scales with how
much of the target the design covers and how well the two match, which means
truncating or mismatching costs binding in roughly the right direction. **The
magnitudes are made up and should not be used for a design decision.** It exists
so the pipeline, the schema and the queries can be exercised in CI without a
licence.

## Layout

```
config/         sweep grid and ranking thresholds
data/           the scan CSV and the generated database
sql/            generated DDL, the two ranking queries, and the trend template
scripts/        DDL rendering, notebook generation, placeholder catalog
src/
  db.py         engine, session, SQLite pragmas, ANALYZE, DDL rendering
  models.py     the schema, single source of truth
  engines.py    ThermoEngine protocol, NUPACK adapter, stub
  catalog.py    scan CSV to variants table, with derived columns
  importer.py   wide energy columns to long results rows
  sweep.py      grid expansion, pending cells, resumable runner
  selection.py  config to bind parameters, query execution, persistence
  trends.py     whitelisted dimensions, trend query assembly
  cli.py        entry point
notebooks/      report notebook, generated so diffs stay readable
tests/          60 tests against an in-memory database
```

## Tests

```
$ make test
60 passed
```

They run against an in-memory database built from the same models, with a
9 variant miniature of the scan: one 12 nt parent, three truncation
combinations, an unmutated reference and two double mutants at each. Enough for
the per-truncation ddG join and the per-class cap to have something to work on.

The ones worth knowing about:

- ranking is monotonically decreasing, and the mirror query returns the exact
  reverse order
- both query files take an identical parameter list
- ddG uses the reference at the matching truncation, and reference rows come
  back at exactly zero
- `delta_g_binding` equals `complex_dg - cdna_dg - target_dg` for every stored
  row, so the denormalised column cannot drift from its components
- the target strand is evaluated once per condition, not once per variant
- rerunning a finished sweep computes zero cells, and widening the grid computes
  exactly the difference
- every trend dimension executes, group counts match the catalog, and an unknown
  dimension name raises rather than reaching SQL
- the generated DDL matches the models

## Notes

`data/scan_25nt.csv` is the 25 nt cDNA double mutation and truncation scan.
`load-variants` needs only a `sequence` column; the scan metadata columns are
optional and any energy columns are ignored at that stage. Variant names are
derived from the scan coordinates (`t5_4_T5G,G9C`), so reloading a re-sorted
file matches existing rows instead of duplicating them. A file with no scan
metadata gets a short sequence hash appended instead.

`length` and `gc_content` are recomputed from the sequence and cross checked
against the file rather than trusted, so a stale column cannot quietly poison a
trend query.

SQLite is used because this is a single-writer batch job and a file that can be
copied around is convenient. Nothing in the schema is SQLite specific apart from
`INTEGER PRIMARY KEY`. The window function in the ranking query needs SQLite
3.25 or newer, which `db.py` checks at connect time.

One gotcha if you edit any of the SQL: SQLAlchemy's `text()` scans the whole
file for bind markers and does not skip `--` comments, so a parameter name
written with its leading colon inside a comment registers as an extra parameter
and the execute fails on a binding count mismatch. The parameter lists at the
top of those files are written without colons for this reason.
