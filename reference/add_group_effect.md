# Add a Designed Nominal Group Effect to a Dataset

Adds a categorical grouping factor (with named levels and chosen
proportions) to a tabular dataset and shifts numeric features by a
\*\*per-level x per-feature\*\* amount, in SD units. Unlike a
Gaussian-copula phenotype (see \[simulate_dataset()\]'s \`inject\`,
which models an \*ordinal\* categorical along a single latent gradient),
this is a genuinely \*\*nominal\*\* effect applied post-hoc: each level
moves its own subset of features in its own directions, so the factor
cannot be reduced to a single ordering.

## Usage

``` r
add_group_effect(
  x,
  name,
  levels,
  probs = NULL,
  effects = "spike_slab",
  effect_prop = 0.1,
  effect_sd = 0.5,
  on = "all",
  reference_level = NULL,
  group_col = name,
  seed = NULL
)
```

## Arguments

- x:

  A \`data.frame\` or a \`synthetica_sim\` (the return value of
  \[simulate_dataset()\]). For a \`synthetica_sim\`, \`\$data\` is
  modified.

- name:

  Character, a label for the grouping factor (used as the default
  \`group_col\` and recorded in the bookkeeping).

- levels:

  Character vector (length \>= 2, unique) of level names.

- probs:

  Optional numeric vector of length \`levels\` summing to 1 giving the
  level proportions. \`NULL\` (default) means equal-sized groups.

- effects:

  Either the string \`"spike_slab"\` (generative, default) or a named
  list keyed by level, each a named numeric vector of per-feature shifts
  in SD units (explicit).

- effect_prop:

  Generative only: per-feature probability that a level perturbs it.
  Default 0.1.

- effect_sd:

  Generative only: SD of the \`N(0, effect_sd)\` slab (the effect-size
  distribution), in each feature's SD units. Default 0.5.

- on:

  Generative only: which features are eligible – \`"all"\` numeric
  columns, a numeric fraction in \`(0, 1\]\` (a random subset), or a
  character vector of feature names. Default \`"all"\`.

- reference_level:

  Optional level name left unshifted (a clean "reference" group).
  \`NULL\` (default) perturbs every level.

- group_col:

  Character, the name of the new column. Must not collide with existing
  columns. Defaults to \`name\`.

- seed:

  Optional integer for reproducible level assignment and (in generative
  mode) the spike-slab draws.

## Value

The same type as \`x\`. For a \`data.frame\`: the input with the
grouping column appended and target features shifted, carrying an
\`attr(., "group_effect")\` record. For a \`synthetica_sim\`: the same
object with \`\$data\` modified and a \`\$group_effect\` slot recording
the realized shift matrix (levels x targets), the assignment, \`probs\`,
and the model parameters.

## Details

Two ways to specify the effect:

\* \*\*Generative\*\* (\`effects = "spike_slab"\`, the default) – for
each non-reference level, each eligible feature is hit with probability
\`effect_prop\` and, if hit, shifted by \`N(0, effect_sd)\` SD units.
Designed for high-dimensional 'omics matrices (Olink, SomaScan,
Metabolon) where an explicit per-feature matrix is infeasible. \*
\*\*Explicit\*\* – \`effects\` is a named list keyed by level, each
element a named numeric vector of per-feature shifts (SD units). Omitted
level = reference (no shift); omitted feature = 0. Designed for
low-dimensional, hand-tuned effects.

Because the shift changes group means, each affected feature's marginal
becomes a mixture across groups – this is intended (the group genuinely
changes the distribution) and is why the effect is applied post-hoc
rather than folded into the copula.

## See also

\[add_batch_effect()\], \[simulate_dataset()\]

## Examples

``` r
set.seed(1)
df <- as.data.frame(matrix(rnorm(300 * 50), 300, 50,
                           dimnames = list(NULL, paste0("p", 1:50))))

# Generative: a 3-level "tissue" factor; each level perturbs ~15% of the
# 50 features, effects ~ N(0, 0.5 SD). Scales directly to 'omics widths.
g <- add_group_effect(df, name = "tissue",
                      levels = c("liver", "muscle", "fat"),
                      probs  = c(0.5, 0.3, 0.2),
                      effect_prop = 0.15, effect_sd = 0.5, seed = 42)
table(g$tissue)
#> 
#>  liver muscle    fat 
#>    150     90     60 
rowSums(attr(g, "group_effect")$shift_matrix != 0)  # features hit per level
#>  liver muscle    fat 
#>      3      9      4 

# Explicit: hand-tuned effects on a couple of named features.
g2 <- add_group_effect(df, name = "soil",
                       levels = c("clay", "sand", "silt"),
                       effects = list(clay = c(p1 = -0.8, p2 = -0.6),
                                      silt = c(p1 =  0.9)),
                       seed = 1)
```
