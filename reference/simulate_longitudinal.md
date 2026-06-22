# Simulate Longitudinal / Repeated-Measures Data

Generate privacy-preserving synthetic copies of longitudinal (panel)
data – repeated measures on the same subjects over time – by identifying
a subject \`id\` column and a \`time\` column (e.g. \`Chick\` and
\`Time\` in \`datasets::ChickWeight\`).

## Usage

``` r
simulate_longitudinal(
  data,
  id,
  time,
  n = NULL,
  static = NULL,
  id_prefix = "subj",
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A long-format \`data.frame\` with one row per (subject, time).

- id:

  Character, the name of the subject identifier column.

- time:

  Character, the name of the time column. Its unique observed values
  define the time grid; genuinely irregular/continuous times should be
  binned to a common grid before calling.

- n:

  Number of synthetic \*\*subjects\*\* to generate. Default \`NULL\`
  uses the number of input subjects.

- static:

  Optional character vector of time-invariant covariate columns
  (constant within each subject, e.g. \`Diet\`). \`NULL\` (default)
  auto-detects them. Columns declared static that actually vary within a
  subject emit a warning and use the first value per subject.

- id_prefix:

  Character prefix for the synthetic subject ids. Default \`"subj"\`
  (ids become \`subj1\`, \`subj2\`, ...).

- verbose:

  Logical; passed to \[simulate_dataset()\].

- ...:

  Further arguments forwarded to \[simulate_dataset()\] (\`seed\`,
  \`privacy\`, \`missingness\`, \`marginal_resolution\`,
  \`max_ordinal_levels\`, ...).

## Value

An object of class \`"synthetica_long"\`: a list with \* \`data\` – the
synthetic long-format \`data.frame\` (same columns as input) \* \`wide\`
– the underlying \`synthetica_sim\` fitted on the wide frame \* \`id\`,
\`time\`, \`times\`, \`static\`, \`features\`, \`n_subjects\` – metadata

## Details

The data are reshaped to \*\*wide\*\* format (one row per subject, with
a column for each time-varying feature at each timepoint), passed
through \[simulate_dataset()\], and reshaped back to long. This turns
within-subject temporal autocorrelation into ordinary between-column
correlation, so the full copula machinery applies for free:
per-timepoint marginals capture the trajectory shape (e.g. growth), the
missingness model reproduces dropout, time-invariant covariates ride
along as single columns, and integer / discrete / polychoric handling
all work.

This handles ChickWeight-style panels, before/after RCTs (two
timepoints), crossover trials (columns per period), and modest biomarker
panels. It does \*\*not\*\* scale to high-dimensional longitudinal
'omics (a warning fires when \`features x timepoints\` exceeds a few
hundred columns), which needs a different method.

## See also

\[simulate_dataset()\]

## Examples

``` r
data(ChickWeight, package = "datasets")
sim <- simulate_longitudinal(ChickWeight, id = "Chick", time = "Time",
                             seed = 1, verbose = FALSE)
head(sim$data)
#>   weight Time Chick Diet
#> 1     43    0 subj1    1
#> 2     51    2 subj1    1
#> 3     61    4 subj1    1
#> 4     72    6 subj1    1
#> 5     84    8 subj1    1
#> 6     87   10 subj1    1
# synthetic chicks have growing weight trajectories and a Diet covariate
```
