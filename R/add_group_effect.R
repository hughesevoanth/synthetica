#' Add a Designed Nominal Group Effect to a Dataset
#'
#' Adds a categorical grouping factor (with named levels and chosen
#' proportions) to a tabular dataset and shifts numeric features by a
#' **per-level x per-feature** amount, in SD units. Unlike a Gaussian-copula
#' phenotype (see [simulate_dataset()]'s `inject`, which models an *ordinal*
#' categorical along a single latent gradient), this is a genuinely **nominal**
#' effect applied post-hoc: each level moves its own subset of features in its
#' own directions, so the factor cannot be reduced to a single ordering.
#'
#' Two ways to specify the effect:
#'
#' * **Generative** (`effects = "spike_slab"`, the default) -- for each
#'   non-reference level, each eligible feature is hit with probability
#'   `effect_prop` and, if hit, shifted by `N(0, effect_sd)` SD units. Designed
#'   for high-dimensional 'omics matrices (Olink, SomaScan, Metabolon) where an
#'   explicit per-feature matrix is infeasible.
#' * **Explicit** -- `effects` is a named list keyed by level, each element a
#'   named numeric vector of per-feature shifts (SD units). Omitted level =
#'   reference (no shift); omitted feature = 0. Designed for low-dimensional,
#'   hand-tuned effects.
#'
#' Because the shift changes group means, each affected feature's marginal
#' becomes a mixture across groups -- this is intended (the group genuinely
#' changes the distribution) and is why the effect is applied post-hoc rather
#' than folded into the copula.
#'
#' @param x A `data.frame` or a `synthetica_sim` (the return value of
#'   [simulate_dataset()]). For a `synthetica_sim`, `$data` is modified.
#' @param name Character, a label for the grouping factor (used as the default
#'   `group_col` and recorded in the bookkeeping).
#' @param levels Character vector (length >= 2, unique) of level names.
#' @param probs Optional numeric vector of length `levels` summing to 1 giving
#'   the level proportions. `NULL` (default) means equal-sized groups.
#' @param effects Either the string `"spike_slab"` (generative, default) or a
#'   named list keyed by level, each a named numeric vector of per-feature
#'   shifts in SD units (explicit).
#' @param effect_prop Generative only: per-feature probability that a level
#'   perturbs it. Default 0.1.
#' @param effect_sd Generative only: SD of the `N(0, effect_sd)` slab (the
#'   effect-size distribution), in each feature's SD units. Default 0.5.
#' @param on Generative only: which features are eligible -- `"all"` numeric
#'   columns, a numeric fraction in `(0, 1]` (a random subset), or a character
#'   vector of feature names. Default `"all"`.
#' @param reference_level Optional level name left unshifted (a clean
#'   "reference" group). `NULL` (default) perturbs every level.
#' @param group_col Character, the name of the new column. Must not collide with
#'   existing columns. Defaults to `name`.
#' @param seed Optional integer for reproducible level assignment and (in
#'   generative mode) the spike-slab draws.
#'
#' @return The same type as `x`. For a `data.frame`: the input with the grouping
#'   column appended and target features shifted, carrying an
#'   `attr(., "group_effect")` record. For a `synthetica_sim`: the same object
#'   with `$data` modified and a `$group_effect` slot recording the realized
#'   shift matrix (levels x targets), the assignment, `probs`, and the model
#'   parameters.
#'
#' @examples
#' set.seed(1)
#' df <- as.data.frame(matrix(rnorm(300 * 50), 300, 50,
#'                            dimnames = list(NULL, paste0("p", 1:50))))
#'
#' # Generative: a 3-level "tissue" factor; each level perturbs ~15% of the
#' # 50 features, effects ~ N(0, 0.5 SD). Scales directly to 'omics widths.
#' g <- add_group_effect(df, name = "tissue",
#'                       levels = c("liver", "muscle", "fat"),
#'                       probs  = c(0.5, 0.3, 0.2),
#'                       effect_prop = 0.15, effect_sd = 0.5, seed = 42)
#' table(g$tissue)
#' rowSums(attr(g, "group_effect")$shift_matrix != 0)  # features hit per level
#'
#' # Explicit: hand-tuned effects on a couple of named features.
#' g2 <- add_group_effect(df, name = "soil",
#'                        levels = c("clay", "sand", "silt"),
#'                        effects = list(clay = c(p1 = -0.8, p2 = -0.6),
#'                                       silt = c(p1 =  0.9)),
#'                        seed = 1)
#'
#' @seealso [add_batch_effect()], [simulate_dataset()]
#' @export
#' @importFrom stats sd
add_group_effect <- function(x, name, levels, probs = NULL,
                             effects         = "spike_slab",
                             effect_prop     = 0.1,
                             effect_sd       = 0.5,
                             on              = "all",
                             reference_level = NULL,
                             group_col       = name,
                             seed            = NULL) {
  is_sim <- inherits(x, "synthetica_sim")
  df     <- if (is_sim) x$data else x

  stopifnot("x must be a data.frame or synthetica_sim" = is.data.frame(df))
  stopifnot("name must be a single string" =
              is.character(name) && length(name) == 1L)
  stopifnot("levels must be a character vector of length >= 2" =
              is.character(levels) && length(levels) >= 2L)
  stopifnot("levels must be unique" = !anyDuplicated(levels))
  k <- length(levels)

  if (is.null(probs)) probs <- rep(1 / k, k)
  stopifnot("probs must be numeric of length(levels), all > 0, summing to 1" =
              is.numeric(probs) && length(probs) == k && all(probs > 0) &&
              abs(sum(probs) - 1) < 1e-6)

  if (group_col %in% names(df)) {
    stop("group_col '", group_col, "' collides with an existing column")
  }

  ref_idx <- NULL
  if (!is.null(reference_level)) {
    stopifnot("reference_level must be one of `levels`" = reference_level %in% levels)
    ref_idx <- match(reference_level, levels)
  }

  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  generative   <- is.character(effects) && length(effects) == 1L &&
    identical(effects, "spike_slab")

  if (generative) {
    stopifnot("effect_prop must be a single number in [0, 1]" =
                is.numeric(effect_prop) && length(effect_prop) == 1L &&
                effect_prop >= 0 && effect_prop <= 1)
    stopifnot("effect_sd must be a single positive number" =
                is.numeric(effect_sd) && length(effect_sd) == 1L && effect_sd > 0)
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
    shift_matrix <- .build_shift_matrix(TRUE, NULL, targets, k, levels,
                                        effect_prop, effect_sd, ref_idx, seed)
  } else {
    stopifnot("effects must be \"spike_slab\" or a named list" =
                is.list(effects) && length(effects) >= 1L && !is.null(names(effects)))
    bad_lv <- setdiff(names(effects), levels)
    if (length(bad_lv)) {
      stop("effects names not in `levels`: ", paste(bad_lv, collapse = ", "))
    }
    targets <- unique(unlist(lapply(effects, names)))
    bad_f   <- setdiff(targets, numeric_cols)
    if (length(bad_f)) {
      stop("effects reference non-numeric or absent feature(s): ",
           paste(bad_f, collapse = ", "))
    }
    shift_matrix <- matrix(0, nrow = k, ncol = length(targets),
                           dimnames = list(levels, targets))
    for (lv in names(effects)) {
      e <- effects[[lv]]
      stopifnot("each effects entry must be a named numeric vector" =
                  is.numeric(e) && !is.null(names(e)))
      for (f in names(e)) shift_matrix[lv, f] <- e[[f]]
    }
  }

  if (length(targets) == 0L) {
    warning("no numeric target features to apply group effect to")
  }

  # Assign rows to levels by proportion (shuffled), then shift in SD units.
  if (!is.null(seed)) set.seed(seed + 1L)
  n      <- nrow(df)
  sizes  <- .batch_sizes(n, k, probs)
  assign <- sample(rep(seq_len(k), times = sizes))

  for (j in targets) {
    s <- stats::sd(df[[j]], na.rm = TRUE)
    if (is.na(s) || s == 0) next
    was_int <- is.integer(df[[j]])
    col     <- df[[j]]
    for (l in seq_len(k)) {
      idx <- assign == l
      col[idx] <- col[idx] + shift_matrix[l, j] * s
    }
    df[[j]] <- if (was_int) .to_int(col) else col
  }

  df[[group_col]] <- factor(levels[assign], levels = levels)

  effect <- list(
    name            = name,
    levels          = levels,
    probs           = probs,
    mode            = if (generative) "spike_slab" else "explicit",
    shift_matrix    = shift_matrix,
    targets         = targets,
    reference_level = reference_level,
    effect_prop     = if (generative) effect_prop else NULL,
    effect_sd       = if (generative) effect_sd else NULL,
    assignment      = assign
  )

  if (is_sim) {
    x$data         <- df
    x$group_effect <- effect
    return(x)
  }
  attr(df, "group_effect") <- effect
  df
}
