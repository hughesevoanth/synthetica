#' Add a Synthetic Batch Effect to a Dataset
#'
#' Adds batch labels to a tabular dataset and (optionally) applies per-batch
#' shifts -- in SD units of each target feature -- to simulate technical batch
#' effects that are common in metabolomics, proteomics, and other 'omics
#' pipelines.
#'
#' Two shift models are supported. A **uniform** model (pass `shifts` as a
#' numeric vector of length `n_batches`) moves every target feature in a batch
#' by the same SD-multiple. A **spike-and-slab** model (pass
#' `shifts = "spike_slab"`) draws a per-feature shift for each batch so that
#' most features are unaffected and a minority are strongly perturbed -- the
#' realistic, heterogeneous signature of a real batch effect.
#'
#' This is a small post-processing helper that operates on either a plain
#' `data.frame` or the output of [simulate_dataset()]. It's intended for
#' generating teaching/test data with a known batch structure and for
#' stress-testing batch-correction methods.
#'
#' @param x A `data.frame` or a `synthetica_sim` object (the return value of
#'   [simulate_dataset()]). When a `synthetica_sim` is passed, the modification
#'   is applied to its `$data` slot.
#' @param n_batches Integer, the number of batches. Default 4.
#' @param shifts Controls the per-batch mean shift, in **SD units** of each
#'   target feature. One of:
#'   * `NULL` (default) -- all zeros; a batch label is added but feature values
#'     are unchanged.
#'   * a numeric vector of length `n_batches` -- the *uniform* model: every
#'     target feature in batch `b` is shifted by `shifts[b] * sd(feature)`.
#'   * the string `"spike_slab"` -- the *generative* model: for each
#'     non-reference batch and each target feature, a shift is drawn as
#'     `rbinom(1, 1, shift_prob) * rnorm(1, 0, shift_sd)`. See `shift_prob` and
#'     `shift_sd`.
#' @param shift_prob Numeric in `[0, 1]`, the per-feature probability of being
#'   "hit" by a shift under the `"spike_slab"` model (the slab weight).
#'   Default 0.2. Ignored for the uniform model.
#' @param shift_sd Positive numeric, the SD of the slab `N(0, shift_sd)` under
#'   the `"spike_slab"` model, in each feature's SD units. Default 0.5. Ignored
#'   for the uniform model.
#' @param reference_batch Integer index of a batch to leave unshifted as a
#'   clean anchor under the `"spike_slab"` model (the reference you correct
#'   *toward*). Default `1L`. `NULL` perturbs all batches. Ignored for the
#'   uniform model.
#' @param on One of `"all"` (every numeric column), a numeric fraction in
#'   `(0, 1]` (a random subset of numeric columns), or a character vector
#'   of feature names. Categorical and non-numeric columns are ignored.
#' @param proportions Optional numeric vector of length `n_batches` summing to
#'   1, giving the relative size of each batch. `NULL` (default) means equally
#'   sized batches (to within one row).
#' @param batch_col Character, the name of the new column. Must not collide
#'   with existing column names. Default `"batch"`.
#' @param level_names Optional character vector of length `n_batches`. If
#'   `NULL` (default), levels are `"batch1"`, `"batch2"`, ...
#' @param seed Optional integer for reproducible batch assignment, target
#'   selection (when `on` is a fraction), and spike-and-slab shift draws.
#'
#' @details
#' Each row is randomly assigned to a batch. With `proportions = NULL` the
#' batches are equally sized (trailing rows distributed one-per-batch when `n`
#' is not divisible by `n_batches`); otherwise sizes follow `proportions` via
#' largest-remainder rounding so they still sum to `n`. For each target feature
#' `f` and batch `b`, values are shifted by `shift[b, f] * sd(f, na.rm = TRUE)`,
#' where `shift[b, f]` comes from the chosen model. NAs are left as NA.
#' Non-numeric columns are passed through unchanged.
#'
#' @return The same type as `x`. For a `data.frame`: the input with the batch
#'   column appended and target features shifted, carrying an
#'   `attr(., "batch_effect")` record. For a `synthetica_sim`: the same object
#'   with `$data` modified and a `$batch_effect` slot recording the assignment
#'   vector, the realized shift matrix (batches x targets), the targets,
#'   `proportions`, and the model parameters.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   x = rnorm(400, 50, 10),
#'   y = rlnorm(400, log(20), 0.5),
#'   z = rnorm(400, 0, 1)
#' )
#'
#' # Uniform model: three indistinguishable batches and one with a +0.5 SD
#' # shift on every numeric feature.
#' df2 <- add_batch_effect(
#'   df,
#'   n_batches = 4,
#'   shifts    = c(0.00, 0.05, -0.05, 0.50),
#'   seed      = 42
#' )
#' table(df2$batch)
#'
#' # Spike-and-slab model: batch 1 is a clean reference; each other batch
#' # perturbs ~20% of features. Unequal batch sizes via `proportions`.
#' df3 <- add_batch_effect(
#'   df,
#'   n_batches   = 4,
#'   shifts      = "spike_slab",
#'   shift_prob  = 0.2,
#'   shift_sd    = 0.5,
#'   proportions = c(0.4, 0.3, 0.2, 0.1),
#'   seed        = 42
#' )
#' attr(df3, "batch_effect")$shift_matrix
#'
#' @seealso [simulate_dataset()]
#' @export
#' @importFrom stats sd rnorm rbinom
add_batch_effect <- function(x,
                             n_batches       = 4L,
                             shifts          = NULL,
                             shift_prob      = 0.2,
                             shift_sd        = 0.5,
                             reference_batch = 1L,
                             on              = "all",
                             proportions     = NULL,
                             batch_col       = "batch",
                             level_names     = NULL,
                             seed            = NULL) {
  is_sim <- inherits(x, "synthetica_sim")
  df     <- if (is_sim) x$data else x

  stopifnot("x must be a data.frame or synthetica_sim" = is.data.frame(df))
  stopifnot("n_batches must be a positive integer scalar" =
              is.numeric(n_batches) && length(n_batches) == 1L && n_batches >= 1)
  n_batches <- as.integer(n_batches)

  # Determine the shift model from `shifts`
  spike_slab <- is.character(shifts) && length(shifts) == 1L &&
    identical(shifts, "spike_slab")

  if (spike_slab) {
    stopifnot("shift_prob must be a single number in [0, 1]" =
                is.numeric(shift_prob) && length(shift_prob) == 1L &&
                shift_prob >= 0 && shift_prob <= 1)
    stopifnot("shift_sd must be a single positive number" =
                is.numeric(shift_sd) && length(shift_sd) == 1L && shift_sd > 0)
    if (!is.null(reference_batch)) {
      stopifnot("reference_batch must be a single index in 1:n_batches" =
                  is.numeric(reference_batch) && length(reference_batch) == 1L &&
                  reference_batch >= 1 && reference_batch <= n_batches)
      reference_batch <- as.integer(reference_batch)
    }
  } else {
    if (is.null(shifts)) shifts <- rep(0, n_batches)
    stopifnot("shifts must be NULL, \"spike_slab\", or numeric of length n_batches" =
                is.numeric(shifts) && length(shifts) == n_batches)
  }

  if (is.null(level_names)) level_names <- paste0("batch", seq_len(n_batches))
  stopifnot("level_names must be character of length n_batches" =
              is.character(level_names) && length(level_names) == n_batches)

  if (!is.null(proportions)) {
    stopifnot("proportions must be numeric of length n_batches, all > 0" =
                is.numeric(proportions) && length(proportions) == n_batches &&
                all(proportions > 0))
    stopifnot("proportions must sum to 1" =
                abs(sum(proportions) - 1) < 1e-6)
  }

  if (batch_col %in% names(df)) {
    stop("batch_col '", batch_col, "' collides with an existing column")
  }

  # Resolve `on` to numeric target columns
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  targets <- if (identical(on, "all")) {
    numeric_cols
  } else if (is.numeric(on) && length(on) == 1L && on > 0 && on <= 1) {
    n_pick <- max(1L, round(on * length(numeric_cols)))
    if (!is.null(seed)) set.seed(seed)
    sample(numeric_cols, n_pick)
  } else if (is.character(on)) {
    intersect(on, numeric_cols)
  } else {
    stop("`on` must be 'all', a numeric fraction in (0,1], or a character vector")
  }
  if (length(targets) == 0L) {
    warning("no numeric target columns to apply batch effect to")
  }

  # Assign batches
  if (!is.null(seed)) set.seed(seed + 1L)
  n     <- nrow(df)
  sizes <- .batch_sizes(n, n_batches, proportions)
  batch <- sample(rep(seq_len(n_batches), times = sizes))

  # Build the per-batch x per-target shift matrix (in SD units)
  shift_matrix <- .build_shift_matrix(
    spike_slab, shifts, targets, n_batches, level_names,
    shift_prob, shift_sd, reference_batch, seed
  )

  # Apply shifts in SD units of each feature
  for (j in targets) {
    s <- stats::sd(df[[j]], na.rm = TRUE)
    if (is.na(s) || s == 0) next
    for (b in seq_len(n_batches)) {
      idx <- batch == b
      df[idx, j] <- df[idx, j] + shift_matrix[b, j] * s
    }
  }

  df[[batch_col]] <- factor(level_names[batch], levels = level_names)

  effect <- list(
    n_batches       = n_batches,
    model           = if (spike_slab) "spike_slab" else "uniform",
    shifts          = if (spike_slab) NULL else shifts,
    shift_matrix    = shift_matrix,
    targets         = targets,
    level_names     = level_names,
    proportions     = proportions,
    shift_prob      = if (spike_slab) shift_prob else NULL,
    shift_sd        = if (spike_slab) shift_sd else NULL,
    reference_batch = if (spike_slab) reference_batch else NULL,
    assignment      = batch
  )

  if (is_sim) {
    x$data         <- df
    x$batch_effect <- effect
    return(x)
  }
  attr(df, "batch_effect") <- effect
  df
}

# Allocate `n` rows across `n_batches` batches. With proportions = NULL the
# batches are equal (trailing rows spread one-per-batch); otherwise sizes follow
# proportions via largest-remainder rounding so they sum to exactly n.
.batch_sizes <- function(n, n_batches, proportions = NULL) {
  if (is.null(proportions)) {
    base <- rep(n %/% n_batches, n_batches)
    rem  <- n %% n_batches
    if (rem > 0) base[seq_len(rem)] <- base[seq_len(rem)] + 1L
    return(base)
  }
  raw   <- proportions * n
  floored <- floor(raw)
  rem   <- n - sum(floored)
  if (rem > 0) {
    # hand the leftover rows to the largest fractional parts
    order_frac <- order(raw - floored, decreasing = TRUE)
    floored[order_frac[seq_len(rem)]] <- floored[order_frac[seq_len(rem)]] + 1L
  }
  as.integer(floored)
}

# Construct the n_batches x n_targets shift matrix (SD units), row/col named.
.build_shift_matrix <- function(spike_slab, shifts, targets, n_batches,
                                level_names, shift_prob, shift_sd,
                                reference_batch, seed) {
  nt <- length(targets)
  m  <- matrix(0, nrow = n_batches, ncol = nt,
               dimnames = list(level_names, targets))
  if (nt == 0L) return(m)

  if (!spike_slab) {
    # uniform model: every target gets the per-batch scalar
    m[] <- matrix(rep(shifts, times = nt), nrow = n_batches, ncol = nt)
    return(m)
  }

  # spike-and-slab: draw per-batch, per-feature; reference batch stays at 0
  if (!is.null(seed)) set.seed(seed + 2L)
  for (b in seq_len(n_batches)) {
    if (!is.null(reference_batch) && b == reference_batch) next
    hit <- stats::rbinom(nt, size = 1L, prob = shift_prob)
    m[b, ] <- hit * stats::rnorm(nt, mean = 0, sd = shift_sd)
  }
  m
}
