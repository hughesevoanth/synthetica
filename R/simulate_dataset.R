# =============================================================================
# synthetica -- simulate_dataset()
#
# Gaussian-copula simulator for mixed-type tabular data.
#
#   Baseline simulation: synthetic copy of a real dataset, preserving
#     marginal distributions, missingness, and joint correlations.
#   Phenotype injection: append synthetic phenotype columns with a
#     user-specified marginal and effect size on a chosen subset of features.
#
# Known limitations:
#   1. Binary phenotypes retain a single thresholding factor on the round-trip:
#      the generating latent is spliced into the copula so the designed
#      correlation is captured exactly, leaving only phi(c) / sqrt(p*(1-p)),
#      where c = qnorm(1-p) and p is the binary prevalence. This factor is
#      reported per binary phenotype in inject_truth$<name>$binary_sim_factor;
#      expected_sim_cor = latent_target_cor * binary_sim_factor.
#   2. Quantitative phenotypes that target many highly-correlated features
#      will see mild "spread-thin" attenuation because R_X^{-1} via ginv is
#      rank-deficient when the target set is dense. The realized |cor| stored
#      in inject_truth lets the caller see this directly.
#   3. Unordered categorical features with k > 2 levels introduce ordering-
#      dependent latent correlations. Marginal level frequencies are preserved
#      exactly; cross-correlations involving such columns are approximate.
# =============================================================================

# Clamp for qnorm(u) to avoid +/- Inf when u hits 0 or 1 exactly.
.SIM_EPS <- 1e-6

# null-coalesce
`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# .detect_col_types
# Return a character vector of column types: "quantitative" / "binary" / "categorical".
# -----------------------------------------------------------------------------
.detect_col_types <- function(df, override = NULL) {
  types <- vapply(names(df), function(nm) {
    x <- df[[nm]]
    if (is.logical(x)) return("binary")
    if (is.factor(x) || is.character(x)) return("categorical")
    if (is.numeric(x)) {
      u <- unique(stats::na.omit(x))
      if (length(u) == 2L) return("binary")
      if (length(u) >= 1L && all(u == round(u))) return("integer")
      return("quantitative")
    }
    stop(sprintf("column '%s' has unsupported type: %s", nm, class(x)[1L]))
  }, character(1L))
  if (!is.null(override)) {
    bad <- setdiff(names(override), names(df))
    if (length(bad)) stop("col_types names not in data: ", paste(bad, collapse = ", "))
    types[names(override)] <- unlist(override)
  }
  types
}


# -----------------------------------------------------------------------------
# .apply_privacy_guardrails
# Strict-mode guardrails:
#   - drop columns with < 2 distinct non-NA values (drop_near_constant)
#   - drop columns with missingness rate > max_missingness
#   - collapse categorical levels with count < min_cell_size into "other"
#   - if min_cell_size < 5, require force = TRUE (loud warning)
# -----------------------------------------------------------------------------
.apply_privacy_guardrails <- function(df, types,
                                      min_cell_size,
                                      drop_near_constant,
                                      max_missingness,
                                      force,
                                      verbose) {
  if (min_cell_size < 5 && !isTRUE(force)) {
    stop("min_cell_size < 5 is below the safe default; pass force = TRUE to override.")
  }
  if (min_cell_size < 5) {
    warning("min_cell_size < 5 -- small categorical cells may permit re-identification.")
  }
  stopifnot("max_missingness must be in [0, 1]" =
              is.numeric(max_missingness) && length(max_missingness) == 1L &&
              max_missingness >= 0 && max_missingness <= 1)

  dropped_const <- character()
  dropped_miss  <- character()
  collapsed     <- list()

  keep <- vapply(names(df), function(nm) {
    x <- df[[nm]]
    n_miss <- mean(is.na(x))
    n_uniq <- length(unique(stats::na.omit(x)))
    const_fail <- isTRUE(drop_near_constant) && n_uniq < 2L
    miss_fail  <- n_miss > max_missingness
    if (const_fail) dropped_const[length(dropped_const) + 1L] <<- nm
    if (miss_fail)  dropped_miss[length(dropped_miss) + 1L]   <<- nm
    !(const_fail || miss_fail)
  }, logical(1L))
  df    <- df[, keep, drop = FALSE]
  types <- types[keep]

  for (nm in names(df)[types == "categorical"]) {
    x   <- as.character(df[[nm]])
    tab <- table(x, useNA = "no")
    rare <- names(tab)[tab < min_cell_size]
    if (length(rare)) {
      x[x %in% rare] <- "other"
      df[[nm]] <- factor(x)
      collapsed[[nm]] <- rare
    } else {
      df[[nm]] <- factor(x)
    }
  }

  audit <- sprintf(
    "privacy: dropped %d near-constant column(s); dropped %d column(s) with missingness > %.0f%%; collapsed rare levels in %d categorical column(s)",
    length(dropped_const), length(dropped_miss), 100 * max_missingness, length(collapsed)
  )
  if (isTRUE(verbose)) message(audit)

  list(data = df, types = types, audit = audit,
       dropped_cols           = c(dropped_const, dropped_miss),
       dropped_near_constant  = dropped_const,
       dropped_high_missing   = dropped_miss,
       collapsed_levels       = collapsed)
}


# -----------------------------------------------------------------------------
# .fit_marginal
# Per-column marginal storage. For quantitative columns we store a quantile
# grid rather than the raw vector -- this is the privacy guarantee.
# -----------------------------------------------------------------------------
# Round to whole numbers and cast to integer when safely in range; otherwise
# keep a rounded double (avoids as.integer() overflow to NA for huge counts).
.to_int <- function(v) {
  v <- round(v)
  mx <- suppressWarnings(max(abs(v), na.rm = TRUE))
  if (!is.finite(mx) || mx <= .Machine$integer.max) as.integer(v) else v
}

.fit_marginal <- function(x, type, resolution,
                          min_cell_size = 10L, max_ordinal_levels = 10L) {
  n     <- length(x)
  obs   <- x[!is.na(x)]
  n_obs <- length(obs)
  if (n_obs < 2L) stop("cannot fit marginal: <2 non-NA observations")

  if (type == "integer") {
    iv     <- as.numeric(obs)
    lv     <- sort(unique(iv))
    counts <- as.numeric(table(factor(iv, levels = lv)))
    ordinal <- length(lv) <= max_ordinal_levels && min(counts) >= min_cell_size
    if (ordinal) {
      return(list(type = "integer", mode = "ordinal",
                  n = n, n_obs = n_obs,
                  levels = lv, probs = counts / sum(counts)))
    }
    probs <- seq(0, 1, length.out = resolution + 1L)
    qs    <- stats::quantile(iv, probs = probs, names = FALSE, type = 7L)
    return(list(type = "integer", mode = "grid",
                n = n, n_obs = n_obs,
                prob_grid = probs, quant_grid = qs,
                mean = mean(iv), sd = stats::sd(iv)))
  }

  if (type == "quantitative") {
    probs <- seq(0, 1, length.out = resolution + 1L)
    qs    <- stats::quantile(obs, probs = probs, names = FALSE, type = 7L)
    return(list(type = "quantitative",
                n = n, n_obs = n_obs,
                prob_grid  = probs,
                quant_grid = qs,
                mean = mean(obs), sd = stats::sd(obs)))
  }

  if (type == "binary") {
    lv <- sort(unique(obs))
    if (length(lv) != 2L) stop("binary column has != 2 distinct values")
    p1 <- mean(obs == lv[2L])
    return(list(type = "binary",
                n = n, n_obs = n_obs,
                levels = lv, probs = c(1 - p1, p1)))
  }

  if (type == "categorical") {
    fx  <- if (is.factor(x)) x else factor(x)
    lv  <- levels(droplevels(fx[!is.na(fx)]))
    tab <- table(factor(obs, levels = lv))
    probs <- as.numeric(tab) / sum(tab)
    return(list(type = "categorical",
                n = n, n_obs = n_obs,
                levels = lv, probs = probs,
                ordered = is.ordered(fx)))
  }

  stop("unknown type: ", type)
}


# -----------------------------------------------------------------------------
# .pit_to_latent
# Probability-integral transform one column to standard normal.
# -----------------------------------------------------------------------------
.pit_to_latent <- function(x, marginal) {
  z <- rep(NA_real_, length(x))
  obs_idx <- which(!is.na(x))
  if (!length(obs_idx)) return(z)

  grid_based <- marginal$type == "quantitative" ||
    (marginal$type == "integer" && marginal$mode == "grid")
  if (grid_based) {
    # ties = mean: silence the "collapsing to unique 'x' values" warning that
    # fires when the quantile grid has duplicates (common when n_obs is small
    # relative to marginal_resolution, or the input has many tied values).
    # The default tie-handling is already mean; passing it explicitly only
    # suppresses the warning.
    u <- stats::approx(marginal$quant_grid, marginal$prob_grid,
                       xout = x[obs_idx], rule = 2L, ties = mean)$y
    u <- pmin(pmax(u, .SIM_EPS), 1 - .SIM_EPS)
    z[obs_idx] <- stats::qnorm(u)
    return(z)
  }

  lv       <- marginal$levels
  probs    <- marginal$probs
  cumprobs <- c(0, cumsum(probs))
  cumprobs[length(cumprobs)] <- 1
  lvl_idx <- match(as.character(x[obs_idx]), as.character(lv))
  if (anyNA(lvl_idx)) stop(".pit_to_latent: value not in stored marginal levels")
  a <- cumprobs[lvl_idx]
  b <- cumprobs[lvl_idx + 1L]
  u <- stats::runif(length(obs_idx), a, b)
  u <- pmin(pmax(u, .SIM_EPS), 1 - .SIM_EPS)
  z[obs_idx] <- stats::qnorm(u)
  z
}


# -----------------------------------------------------------------------------
# .back_transform
# Inverse PIT: latent z to original scale + type via the stored marginal.
# -----------------------------------------------------------------------------
.back_transform <- function(z, marginal) {
  u <- stats::pnorm(z)
  grid_based <- marginal$type == "quantitative" ||
    (marginal$type == "integer" && marginal$mode == "grid")
  if (grid_based) {
    out <- stats::approx(marginal$prob_grid, marginal$quant_grid,
                         xout = u, rule = 2L)$y
    if (marginal$type == "integer") return(.to_int(out))
    return(out)
  }
  lv       <- marginal$levels
  cumprobs <- c(0, cumsum(marginal$probs))
  cumprobs[length(cumprobs)] <- 1
  bin      <- findInterval(u, cumprobs, rightmost.closed = TRUE, all.inside = TRUE)
  out_chr  <- lv[bin]
  if (marginal$type == "integer") return(.to_int(as.numeric(out_chr)))
  if (marginal$type == "binary" && is.numeric(lv)) {
    out <- as.numeric(out_chr)
    if (all(lv == round(lv))) return(.to_int(out))   # whole-valued binary -> integer
    return(out)
  }
  if (marginal$type == "categorical" && isTRUE(marginal$ordered)) {
    return(factor(out_chr, levels = lv, ordered = TRUE))
  }
  if (marginal$type == "categorical") return(factor(out_chr, levels = lv))
  out_chr
}


# -----------------------------------------------------------------------------
# Thresholded-normal latent-correlation recovery for discrete columns.
#
# A discrete column is treated as a latent N(0,1) thresholded at cutpoints.
# These helpers recover the underlying Gaussian correlation that, after
# thresholding, reproduces the observed association -- removing the attenuation
# that the uniform-within-bin PIT would otherwise leave in the copula. Base R
# only (no tetrachoric/polychoric package dependency): binary uses Digby
# (tetrachoric) / biserial; ordered multi-level uses polychoric / polyserial.
# -----------------------------------------------------------------------------

.SIM_CORR_CLIP <- 0.999

# Gauss-Legendre nodes/weights on [-1, 1] via Golub-Welsch (eigen of the
# Legendre Jacobi matrix). Computed once and reused by .bvn_cdf.
.gauss_legendre <- function(n) {
  i  <- seq_len(n - 1L)
  b  <- i / sqrt(4 * i^2 - 1)
  J  <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  ev <- eigen(J, symmetric = TRUE)
  o  <- order(ev$values)
  list(nodes = ev$values[o], weights = 2 * (ev$vectors[1L, o])^2)
}
.GL_NODES <- .gauss_legendre(30L)

# Standard bivariate-normal CDF P(X <= h, Y <= k), corr rho, via the integral
# form Phi(h)Phi(k) + (1/2pi) int_0^rho 1/sqrt(1-t^2) exp(...) dt. Vectorized
# over (h, k); infinite thresholds fall out of the Phi(h)Phi(k) base term.
.bvn_cdf <- function(h, k, rho) {
  base <- stats::pnorm(h) * stats::pnorm(k)
  if (rho == 0) return(base)
  rho <- max(min(rho, 0.999999), -0.999999)
  t   <- rho / 2 * (.GL_NODES$nodes + 1)
  w   <- rho / 2 * .GL_NODES$weights
  omt <- 1 - t^2
  add <- mapply(function(hi, ki) {
    if (!is.finite(hi) || !is.finite(ki)) return(0)
    sum(w / sqrt(omt) * exp(-(hi^2 - 2 * t * hi * ki + ki^2) / (2 * omt)))
  }, h, k)
  base + add / (2 * pi)
}

# Digby's tetrachoric approximation for a pair of 0/1 vectors (0.5 continuity
# correction so empty cells do not give an infinite odds ratio).
.tetrachoric_digby <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(0)
  n11 <- sum(x == 1 & y == 1) + 0.5
  n00 <- sum(x == 0 & y == 0) + 0.5
  n10 <- sum(x == 1 & y == 0) + 0.5
  n01 <- sum(x == 0 & y == 1) + 0.5
  omega <- (n11 * n00) / (n10 * n01)
  rho   <- (omega^(pi / 4) - 1) / (omega^(pi / 4) + 1)
  max(min(rho, .SIM_CORR_CLIP), -.SIM_CORR_CLIP)
}

# Polyserial (moment) estimate between an ordinal code vector (1..k) and a
# continuous vector: rho = r * sd(codes) / sum(phi(thresholds)). Reduces to the
# biserial estimate at k = 2.
.polyserial_latent <- function(codes, cont) {
  ok <- !is.na(codes) & !is.na(cont)
  codes <- codes[ok]; cont <- cont[ok]
  if (length(codes) < 2L) return(0)
  lv <- sort(unique(codes)); k <- length(lv)
  if (k < 2L) return(0)
  p   <- as.numeric(table(factor(codes, levels = lv))) / length(codes)
  tau <- stats::qnorm(cumsum(p)[-k])
  s   <- stats::sd(codes)
  if (is.na(s) || s == 0) return(0)
  r <- suppressWarnings(stats::cor(codes, cont))
  if (is.na(r)) return(0)
  denom <- sum(stats::dnorm(tau))
  if (denom == 0) return(0)
  rho <- r * s / denom
  max(min(rho, .SIM_CORR_CLIP), -.SIM_CORR_CLIP)
}

# Two-step polychoric ML between two ordinal code vectors: thresholds fixed from
# the marginals, then optimise rho over the contingency-table likelihood.
# Generalises tetrachoric (k = m = 2).
.polychoric_latent <- function(ca, cb) {
  ok <- !is.na(ca) & !is.na(cb)
  ca <- ca[ok]; cb <- cb[ok]
  la <- sort(unique(ca)); lb <- sort(unique(cb))
  ka <- length(la); kb <- length(lb)
  if (ka < 2L || kb < 2L) return(0)
  tab <- table(factor(ca, levels = la), factor(cb, levels = lb))
  nij <- matrix(as.numeric(tab), ka, kb)
  n   <- sum(nij)
  ta  <- c(-Inf, stats::qnorm(cumsum(rowSums(nij) / n)[-ka]), Inf)
  tb  <- c(-Inf, stats::qnorm(cumsum(colSums(nij) / n)[-kb]), Inf)
  negll <- function(rho) {
    G <- outer(ta, tb, .bvn_cdf, rho)
    P <- G[-1, -1, drop = FALSE] - G[-(ka + 1L), -1, drop = FALSE] -
         G[-1, -(kb + 1L), drop = FALSE] + G[-(ka + 1L), -(kb + 1L), drop = FALSE]
    -sum(nij * log(pmax(P, 1e-12)))
  }
  rho <- stats::optimize(negll, interval = c(-0.999, 0.999))$minimum
  max(min(rho, .SIM_CORR_CLIP), -.SIM_CORR_CLIP)
}

# Code a discrete column to 0/1 (binary) with 1 = the higher-latent level,
# matching how .fit_marginal / .pit_to_latent order levels.
.binary_code <- function(col) {
  z <- if (is.factor(col)) as.integer(col) else suppressWarnings(as.numeric(col))
  as.integer(z == max(z, na.rm = TRUE))
}

# Classify a column for the discrete correction and return its rank codes.
#   binary2    : exactly 2 observed levels (binary type, or 2-level factor)
#   ordinal    : ordered, > 2 levels (integer ordinal-mode, or ordered factor)
#   continuous : quantitative, or integer grid-mode (faithful via the grid PIT)
#   other      : unordered factor with > 2 levels (order-ambiguous; skipped)
.col_disc_info <- function(col, type, marginal) {
  if (type == "binary") return(list(kind = "binary2", codes = .binary_code(col)))
  if (type == "quantitative") return(list(kind = "continuous", codes = NULL))
  if (type == "integer") {
    if (!is.null(marginal) && identical(marginal$mode, "ordinal")) {
      lv <- sort(unique(stats::na.omit(as.numeric(col))))
      if (length(lv) == 2L) return(list(kind = "binary2", codes = .binary_code(col)))
      return(list(kind = "ordinal", codes = match(as.numeric(col), lv)))
    }
    return(list(kind = "continuous", codes = NULL))   # grid mode: continuous-like
  }
  if (type == "categorical") {
    lv <- if (is.factor(col)) levels(droplevels(col[!is.na(col)]))
          else unique(stats::na.omit(col))
    if (length(lv) == 2L) return(list(kind = "binary2", codes = .binary_code(col)))
    if (is.ordered(col))  return(list(kind = "ordinal",
                                      codes = match(as.character(col), levels(col))))
    return(list(kind = "other", codes = NULL))         # unordered > 2: skipped
  }
  list(kind = "other", codes = NULL)
}

# Replace, in R, the latent-correlation entries for every pair where at least
# one column is an *input* ordered-discrete column (binary or ordinal, origin
# not in `injected`). Injected columns carry their true latent via splicing.
.correct_discrete_corr <- function(R, data, types, marginals, injected = character()) {
  cols <- colnames(R)
  info <- lapply(cols, function(c) .col_disc_info(data[[c]], types[[c]], marginals[[c]]))
  names(info) <- cols
  kind <- vapply(info, `[[`, character(1L), "kind")
  prim <- cols[kind %in% c("binary2", "ordinal") & !(cols %in% injected)]
  if (!length(prim)) return(R)
  changed <- FALSE
  for (a in prim) {
    ia <- info[[a]]
    for (b in cols) {
      if (b == a) next
      ib  <- info[[b]]
      rho <- NULL
      if (ia$kind == "binary2") {
        if (ib$kind == "continuous")   rho <- .polyserial_latent(ia$codes, data[[b]])
        else if (ib$kind == "binary2") rho <- .tetrachoric_digby(ia$codes, ib$codes)
        else if (ib$kind == "ordinal") rho <- .polychoric_latent(ia$codes, ib$codes)
      } else {  # ordinal
        if (ib$kind == "continuous")   rho <- .polyserial_latent(ia$codes, data[[b]])
        else if (ib$kind %in% c("binary2", "ordinal"))
                                       rho <- .polychoric_latent(ia$codes, ib$codes)
      }
      if (!is.null(rho)) {
        R[a, b] <- rho
        R[b, a] <- rho
        changed <- TRUE
      }
    }
  }
  if (changed) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE,
                                  conv.tol = 1e-7, maxit = 200L)$mat)
  }
  R
}

# -----------------------------------------------------------------------------
# .estimate_latent_corr
# Pairwise-complete Pearson correlation on the latent z's, with optional
# Gaussian noise on the off-diagonals, then projected to the nearest PD
# correlation matrix via Matrix::nearPD.
# -----------------------------------------------------------------------------
.estimate_latent_corr <- function(Z, noise_on_corr = 0, min_pair_n = 5L) {
  p <- ncol(Z)
  R <- stats::cor(Z, use = "pairwise.complete.obs")
  M <- !is.na(Z)
  N <- crossprod(M)
  low_overlap <- (N < min_pair_n)
  R[low_overlap | is.na(R)] <- 0
  diag(R) <- 1

  if (noise_on_corr > 0) {
    E <- matrix(stats::rnorm(p * p, sd = noise_on_corr), nrow = p)
    E <- (E + t(E)) / 2
    diag(E) <- 0
    R <- R + E
  }

  R_pd <- as.matrix(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE,
                                   conv.tol = 1e-7, maxit = 200L)$mat)
  R_pd
}


# -----------------------------------------------------------------------------
# .fit_missingness / .sample_missingness
# Independent Bernoulli (MCAR) or Gaussian copula on the NA-indicator matrix.
# -----------------------------------------------------------------------------
.fit_missingness <- function(df, mode) {
  M <- vapply(df, function(x) as.integer(is.na(x)), integer(nrow(df)))
  if (!is.matrix(M)) M <- matrix(M, nrow = nrow(df))
  colnames(M) <- names(df)
  rates <- colMeans(M)

  if (mode == "none") return(list(mode = "none", rates = rep(0, ncol(M)), corr = NULL))
  if (mode == "MCAR") return(list(mode = "MCAR", rates = rates,        corr = NULL))

  active <- which(rates > 0 & rates < 1)
  if (length(active) < 2L) {
    return(list(mode = "preserve_pattern", rates = rates, corr = NULL, active = active))
  }
  marginals_M <- lapply(active, function(j) {
    p <- rates[j]
    list(type = "binary", levels = c(0L, 1L), probs = c(1 - p, p))
  })
  Z <- vapply(seq_along(active), function(k) {
    .pit_to_latent(M[, active[k]], marginals_M[[k]])
  }, numeric(nrow(M)))
  R <- .estimate_latent_corr(Z)
  list(mode = "preserve_pattern", rates = rates, corr = R, active = active)
}


.sample_missingness <- function(fit, n, p, col_names) {
  mask <- matrix(FALSE, nrow = n, ncol = p, dimnames = list(NULL, col_names))
  if (fit$mode == "none") return(mask)

  if (fit$mode == "MCAR" || is.null(fit$corr)) {
    for (j in seq_len(p)) {
      r <- fit$rates[j]
      if (r > 0) mask[, j] <- stats::runif(n) < r
    }
    return(mask)
  }

  active <- fit$active
  Z <- MASS::mvrnorm(n, mu = rep(0, length(active)), Sigma = fit$corr)
  for (k in seq_along(active)) {
    j <- active[k]
    r <- fit$rates[j]
    u <- stats::pnorm(Z[, k])
    mask[, j] <- u >= (1 - r)
  }
  for (j in setdiff(seq_len(p), active)) {
    r <- fit$rates[j]
    if (r == 1) mask[, j] <- TRUE
  }
  mask
}


# -----------------------------------------------------------------------------
# Phenotype-injection helpers
# -----------------------------------------------------------------------------
.validate_inject_spec <- function(inject, data_names) {
  if (!is.list(inject) || is.null(names(inject)) || any(names(inject) == "")) {
    stop("`inject` must be a named list of phenotype specs")
  }
  collisions <- intersect(names(inject), data_names)
  if (length(collisions)) {
    stop("phenotype name(s) collide with feature name(s): ",
         paste(collisions, collapse = ", "))
  }
  valid_kinds <- c("mean_shift", "r_squared")
  valid_types <- c("quantitative", "binary", "categorical")
  for (nm in names(inject)) {
    s <- inject[[nm]]
    if (!is.list(s))                    stop("inject[['", nm, "']] must be a list")
    if (!isTRUE(s$type %in% valid_types))
      stop("inject[['", nm, "']]$type must be one of: ", paste(valid_types, collapse = ", "))
    if (s$type == "quantitative") {
      if (!is.numeric(s$mean) || !is.numeric(s$sd) || s$sd <= 0)
        stop("inject[['", nm, "']] needs numeric mean and positive sd")
    } else if (s$type == "binary") {
      if (!is.numeric(s$prob) || s$prob <= 0 || s$prob >= 1)
        stop("inject[['", nm, "']]$prob must be in (0, 1)")
    } else {  # categorical
      if (!is.character(s$levels) || length(s$levels) < 2L)
        stop("inject[['", nm, "']]$levels must be a character vector of length >= 2")
      if (anyDuplicated(s$levels))
        stop("inject[['", nm, "']]$levels must be unique")
      if (!is.numeric(s$probs) || length(s$probs) != length(s$levels))
        stop("inject[['", nm, "']]$probs must be numeric of the same length as levels")
      if (any(s$probs <= 0) || abs(sum(s$probs) - 1) > 1e-6)
        stop("inject[['", nm, "']]$probs must be positive and sum to 1")
    }
    if (!is.list(s$effect) || !isTRUE(s$effect$kind %in% valid_kinds))
      stop("inject[['", nm, "']]$effect$kind must be one of: ",
           paste(valid_kinds, collapse = ", "))
    if (s$type == "categorical" && s$effect$kind != "r_squared")
      stop("inject[['", nm, "']]: categorical phenotypes support only ",
           "effect$kind = \"r_squared\" in v1")
    if (s$effect$kind == "mean_shift" && !is.numeric(s$effect$size_sd))
      stop("inject[['", nm, "']]$effect$size_sd must be numeric")
    if (s$effect$kind == "r_squared" &&
        (!is.numeric(s$effect$min_r2) || s$effect$min_r2 < 0 || s$effect$min_r2 > 1))
      stop("inject[['", nm, "']]$effect$min_r2 must be in [0, 1]")
    if (is.null(s$effect$on))
      stop("inject[['", nm, "']]$effect$on must be specified")
  }
  inject
}


.resolve_targets <- function(on, data_names) {
  if (identical(on, "all")) return(data_names)
  if (is.numeric(on) && length(on) == 1L && on > 0 && on <= 1) {
    n_pick <- max(1L, round(on * length(data_names)))
    return(sample(data_names, n_pick))
  }
  if (is.character(on)) {
    missing_t <- setdiff(on, data_names)
    if (length(missing_t)) {
      warning("inject target(s) not in data (dropped): ",
              paste(missing_t, collapse = ", "))
    }
    return(intersect(on, data_names))
  }
  stop("`on` must be 'all', a numeric fraction in (0,1], or a character vector")
}


.effect_to_target_corr <- function(phen_spec, n_targets) {
  eff <- phen_spec$effect
  if (eff$kind == "mean_shift") {
    base <- abs(eff$size_sd)
    if (phen_spec$type == "binary") {
      p   <- phen_spec$prob
      phi <- stats::dnorm(stats::qnorm(1 - p))
      base <- base * p * (1 - p) / phi
    }
  } else {  # r_squared
    base <- sqrt(eff$min_r2)
  }
  rep(min(base, 0.999), n_targets)
}


.rank_standardize <- function(X) {
  out <- matrix(0, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (!is.numeric(x)) next
    obs <- !is.na(x)
    if (sum(obs) < 2L) next
    r <- rank(x[obs], ties.method = "average") / (sum(obs) + 1L)
    out[obs, j] <- stats::qnorm(r)
  }
  out
}


.generate_phenotype <- function(phen_spec, X_target, target_corr) {
  X_std <- .rank_standardize(X_target)
  R_X   <- stats::cor(X_std)
  R_X[is.na(R_X)] <- 0
  diag(R_X) <- 1

  w <- as.vector(MASS::ginv(R_X) %*% target_corr)
  quad <- sum(target_corr * w)
  if (quad >= 0.99) {
    lambda <- sqrt(0.99 / quad)
    warning(sprintf(
      "phenotype effect infeasible: requested |corr| mean %.3f, clipped to %.3f",
      mean(abs(target_corr)), mean(abs(target_corr * lambda))
    ))
    target_corr <- target_corr * lambda
    w    <- w    * lambda
    quad <- sum(target_corr * w)
  }
  sigma_eps <- sqrt(max(0, 1 - quad))
  z <- as.vector(X_std %*% w) + stats::rnorm(nrow(X_std), 0, sigma_eps)
  z <- (z - mean(z)) / stats::sd(z)

  if (phen_spec$type == "quantitative") {
    y <- phen_spec$mean + phen_spec$sd * z
  } else if (phen_spec$type == "binary") {
    y <- as.integer(z > stats::qnorm(1 - phen_spec$prob))
  } else {  # categorical: cut the latent at cumulative-probability thresholds
    cutpoints <- stats::qnorm(cumsum(phen_spec$probs))
    cutpoints <- cutpoints[-length(cutpoints)]        # k-1 interior cuts
    bin <- findInterval(z, cutpoints) + 1L            # 1..k; low z -> first level
    y   <- factor(phen_spec$levels[bin], levels = phen_spec$levels)
  }
  # z is the standardized generating latent. It is spliced into the copula's
  # correlation estimation (see simulate_dataset) so the designed latent
  # correlation is captured exactly rather than re-imputed (and attenuated)
  # from the thresholded values.
  list(y = y, w = w, target_corr = target_corr, z = z)
}


.attach_phenotypes <- function(data, types, inject, verbose) {
  inject <- .validate_inject_spec(inject, names(data))
  inject_truth <- list()
  latents      <- list()

  for (nm in names(inject)) {
    s <- inject[[nm]]
    targets <- .resolve_targets(s$effect$on, names(data))
    quant_t <- targets[types[targets] %in% c("quantitative", "binary")]
    if (length(quant_t) < length(targets)) {
      warning("phenotype '", nm, "': dropped ",
              length(targets) - length(quant_t),
              " categorical target(s) -- not supported as effect targets in v1")
    }
    targets <- quant_t
    if (length(targets) == 0L) {
      warning("phenotype '", nm, "': no valid target features remain; skipping")
      next
    }

    target_corr <- .effect_to_target_corr(s, length(targets))
    X_target    <- as.matrix(data[, targets, drop = FALSE])
    gen <- .generate_phenotype(s, X_target, target_corr)

    # categorical phenotypes are factors; correlate their integer (ordinal) code
    y_num <- if (s$type == "categorical") as.integer(gen$y) else gen$y
    real_corr <- vapply(seq_along(targets), function(k) {
      stats::cor(y_num, data[[targets[k]]],
                 use = "pairwise.complete.obs", method = "spearman")
    }, numeric(1L))

    # Single, unavoidable thresholding factor relating the designed latent
    # correlation to the realized point-biserial. The generating latent is
    # spliced into the copula (see simulate_dataset), so the *extra* attenuation
    # that used to square this factor is eliminated.
    bin_factor <- if (s$type == "binary") {
      p   <- s$prob
      phi <- stats::dnorm(stats::qnorm(1 - p))
      phi / sqrt(p * (1 - p))
    } else 1

    latent_target_cor <- mean(abs(gen$target_corr), na.rm = TRUE)

    data[[nm]] <- gen$y
    types[nm]  <- s$type
    latents[[nm]] <- gen$z
    inject_truth[[nm]] <- list(
      spec               = s,
      targets            = targets,
      target_corr        = gen$target_corr,
      realized_cor       = setNames(real_corr, targets),
      mean_abs_real_cor  = mean(abs(real_corr), na.rm = TRUE),
      latent_target_cor  = latent_target_cor,
      binary_sim_factor  = bin_factor,
      expected_sim_cor   = latent_target_cor * bin_factor
    )
    if (isTRUE(verbose)) {
      msg <- sprintf(
        "injected '%s' (%s) onto %d feature(s): realized |Spearman cor| mean = %.3f",
        nm, s$type, length(targets), mean(abs(real_corr), na.rm = TRUE)
      )
      if (s$type == "binary") {
        msg <- paste0(msg, sprintf(
          " (binary: expect sim ~%.3f, thresholding factor %.2f)",
          latent_target_cor * bin_factor, bin_factor
        ))
      }
      message(msg)
    }
  }
  list(data = data, types = types, inject_truth = inject_truth, latents = latents)
}


# -----------------------------------------------------------------------------
# .simulate_stratified
# Per-stratum copula fitting and sampling, driven by simulate_dataset() when
# `stratify_by` is set. Splits the input by the named column, drops it from
# each subset, recursively calls simulate_dataset() (with stratify_by = NULL)
# on each stratum, and concatenates the per-stratum sims with the stratify
# column re-attached. Strata with fewer than `min_strata_n` rows are skipped
# with a warning. Rows of the output stay grouped by stratum.
# -----------------------------------------------------------------------------
.simulate_stratified <- function(data, n, stratify_by, min_strata_n,
                                 col_types, marginal_resolution,
                                 missingness, copula_regularization,
                                 privacy, seed, verbose,
                                 max_ordinal_levels = 10L) {
  strat_col <- data[[stratify_by]]
  if (anyNA(strat_col)) {
    n_na <- sum(is.na(strat_col))
    warning(sprintf("%d row(s) have NA in stratify column '%s' and will be dropped",
                    n_na, stratify_by))
    data      <- data[!is.na(strat_col), , drop = FALSE]
    strat_col <- data[[stratify_by]]
  }

  freq         <- table(strat_col, useNA = "no")
  small_levels <- names(freq)[freq < min_strata_n]
  if (length(small_levels)) {
    warning(sprintf("strata with < %d rows are skipped: %s",
                    min_strata_n, paste(small_levels, collapse = ", ")))
  }
  keep_levels <- setdiff(names(freq), small_levels)
  if (length(keep_levels) < 2L) {
    stop("fewer than 2 strata remain after filtering; cannot stratify")
  }

  # Allocate per-level n proportional to the kept input frequencies; distribute
  # any rounding remainder to the largest strata.
  keep_freq   <- freq[keep_levels]
  n_per_level <- floor(n * as.numeric(keep_freq) / sum(keep_freq))
  remainder   <- n - sum(n_per_level)
  if (remainder > 0L) {
    biggest <- order(as.numeric(keep_freq), decreasing = TRUE)[seq_len(remainder)]
    n_per_level[biggest] <- n_per_level[biggest] + 1L
  }
  names(n_per_level) <- keep_levels

  sim_list       <- list()
  marginals_list <- list()
  corr_list      <- list()
  audit_list     <- list()
  dropped_list   <- list()

  for (k in seq_along(keep_levels)) {
    lev      <- keep_levels[k]
    sub      <- data[strat_col == lev, , drop = FALSE]
    sub[[stratify_by]] <- NULL

    sub_seed <- if (!is.null(seed)) seed + k * 1000L else NULL
    if (isTRUE(verbose)) message(sprintf(
      "stratum '%s': %d input rows -> %d simulated rows", lev, nrow(sub), n_per_level[lev]))

    sub_sim <- simulate_dataset(
      data                  = sub,
      n                     = n_per_level[lev],
      col_types             = col_types,
      marginal_resolution   = marginal_resolution,
      missingness           = missingness,
      inject                = NULL,            # disallowed when stratifying
      copula_regularization = copula_regularization,
      privacy               = privacy,
      seed                  = sub_seed,
      verbose               = FALSE,           # quieter inside the loop
      stratify_by           = NULL,            # prevent infinite recursion
      max_ordinal_levels    = max_ordinal_levels
    )

    sub_sim$data[[stratify_by]] <- lev
    sim_list[[lev]]       <- sub_sim$data
    marginals_list[[lev]] <- sub_sim$marginals
    corr_list[[lev]]      <- sub_sim$corr
    audit_list[[lev]]     <- sub_sim$audit
    dropped_list[[lev]]   <- sub_sim$dropped_cols
  }

  out_data <- do.call(rbind, sim_list)
  rownames(out_data) <- NULL

  if (is.factor(strat_col)) {
    out_data[[stratify_by]] <- factor(out_data[[stratify_by]],
                                      levels = levels(strat_col))
  }
  out_data <- out_data[, names(data), drop = FALSE]

  structure(
    list(
      data                  = out_data,
      marginals             = marginals_list,
      corr                  = corr_list,
      missingness           = NULL,
      inject_truth          = NULL,
      dropped_cols          = unique(unlist(dropped_list)),
      dropped_near_constant = NULL,
      dropped_high_missing  = NULL,
      audit                 = audit_list,
      stratify_by           = stratify_by,
      n_per_level           = n_per_level,
      small_levels_skipped  = small_levels,
      call                  = sys.call(-1L)
    ),
    class = "synthetica_sim"
  )
}


# -----------------------------------------------------------------------------
# simulate_dataset()
# -----------------------------------------------------------------------------

#' Simulate a synthetic dataset preserving marginals and joint correlations
#'
#' Generate a synthetic copy of a mixed-type tabular dataset via a Gaussian
#' copula. The per-column marginals (means, SDs, level frequencies, missingness
#' rate) and the joint correlation structure of the input are preserved
#' without retaining row-level information. Optionally inject synthetic
#' phenotype columns with a user-specified marginal and effect size on a
#' chosen subset of features.
#'
#' @param data A `data.frame` or numeric `matrix`. Samples in rows, features in
#'   columns. Mixed types allowed: numeric (quantitative or binary), logical,
#'   factor, character.
#' @param n Integer, the number of rows in the simulated output. Defaults to
#'   `nrow(data)`.
#' @param col_types Optional named list overriding the auto-detected column
#'   types. Each entry must be one of `"quantitative"`, `"integer"`,
#'   `"binary"`, or `"categorical"`. Whole-number columns are auto-detected as
#'   `"integer"`; override to `"quantitative"` to keep a fractional (float)
#'   representation.
#' @param marginal_resolution Integer, the number of quantile knots used to
#'   represent each quantitative column's marginal CDF. Default 128. Storing
#'   a quantile grid (rather than the raw vector) is the privacy guarantee
#'   for continuous columns.
#' @param missingness One of `"preserve_pattern"` (default), `"MCAR"`, or
#'   `"none"`. With `"preserve_pattern"`, the NA-indicator matrix is modeled
#'   by a second Gaussian copula so that co-missing variables stay co-missing.
#' @param max_ordinal_levels Integer, default 10. An `"integer"` column with at
#'   most this many distinct values (and every value occurring at least
#'   `min_cell_size` times) is modeled as an ordinal trait with exact level
#'   frequencies; otherwise it uses the quantile grid and is rounded to the
#'   nearest integer on output. Either way the output is integer-typed.
#' @param inject Optional named list of phenotype specifications (phenotype injection).
#'   See **Details** and **Examples**.
#' @param copula_regularization Currently unused; reserved for alternate
#'   PD-projection strategies. The default is `Matrix::nearPD()` via
#'   `"nearPD"`.
#' @param privacy A list of privacy knobs:
#'   * `min_cell_size` (default 10): categorical levels with fewer than this
#'     many observations are collapsed into `"other"`. Below 5 requires
#'     `force = TRUE` and emits a warning.
#'   * `drop_near_constant` (default `TRUE`): drop columns with fewer than 2
#'     distinct non-NA values.
#'   * `max_missingness` (default 0.90): drop columns with missingness rate
#'     above this threshold (they cannot be faithfully simulated in their
#'     extreme tail and are a re-identification risk).
#'   * `noise_on_corr` (default 0): standard deviation of Gaussian noise
#'     added to off-diagonals of the latent correlation matrix. Positive
#'     values give a differential-privacy-flavored guarantee; 0 means none.
#'   * `force` (default `FALSE`): set `TRUE` to override the
#'     `min_cell_size < 5` safety check.
#' @param seed Optional integer. When non-`NULL`, seeds the (PIT jitter,
#'   phenotype injection, MVN draw, missingness mask) stages independently
#'   so the output is reproducible. The global RNG is not touched if `NULL`.
#' @param verbose If `TRUE` (default), prints progress and a privacy-audit
#'   line summarising what was dropped or collapsed.
#' @param stratify_by Optional character of length 1 naming a column to
#'   stratify on. When set, the simulator fits a separate copula within each
#'   level of that column and concatenates the per-stratum sims, preserving
#'   within-group joint distributions that a single copula would average
#'   away. Useful when a categorical column drives big mean shifts in
#'   continuous columns (e.g., `Species` in `iris`). Cannot be combined with
#'   `inject` in this version.
#' @param min_strata_n Integer minimum number of rows required per stratum
#'   for fitting. Strata below this are skipped with a warning. Default 30.
#'
#' @details
#' **Method.** The simulator fits a per-column marginal (empirical quantile
#' grid for quantitative columns; level proportions for binary and
#' categorical) and then maps each observation to a standard normal via the
#' probability integral transform (continuity-corrected for discrete
#' variables). Pairwise-complete Pearson correlation of the latent scores
#' is projected to the nearest positive-definite correlation matrix
#' (`Matrix::nearPD()`). Synthetic rows are drawn from `MVN(0, R)` and
#' back-transformed column-wise through the stored marginals. Missingness
#' is modeled by a second Gaussian copula on the NA-indicator matrix.
#'
#' **Phenotype injection.** Each entry of `inject` is a list with
#' fields:
#' * `type` -- `"quantitative"`, `"binary"`, or `"categorical"`
#' * `mean`, `sd` (quantitative); `prob` (binary); **or** `levels` (character,
#'   length >= 2) and `probs` (proportions, same length, all > 0, summing to 1)
#'   for `"categorical"`
#' * `effect = list(kind, ..., on)` where `kind` is `"mean_shift"` (then
#'   `size_sd` gives the standardized effect) or `"r_squared"` (then
#'   `min_r2` gives the target R^2 on the latent scale), and `on` is one of
#'   `"all"`, a numeric fraction in `(0, 1]` selecting a random subset of
#'   target features, or a character vector of feature names. Categorical
#'   phenotypes support only `kind = "r_squared"` in this version.
#'
#' Phenotypes are generated on the rank-standardized latent scale and
#' back-transformed to the user-facing marginal. The feasibility of the
#' requested effect is enforced (`rho' R_X^-1 rho <= 1`); on violation, the
#' targets are scaled down and a warning is emitted.
#'
#' A `"categorical"` phenotype is an **ordinal** trait: a single
#' latent is correlated with the targets and cut at the cumulative-probability
#' thresholds of `probs` into ordered `levels` (the first level is the low end
#' of the gradient). Binary is the two-level special case. Because thresholding
#' attenuates the latent correlation, the realized per-target Spearman
#' association is reported empirically in `inject_truth$<name>$realized_cor`.
#'
#' **Binary fidelity.** Injected binary (and categorical/continuous) phenotypes
#' have their generating latent spliced into the copula, so the designed latent
#' correlation is captured exactly; only the single, unavoidable thresholding
#' factor `phi(c) / sqrt(p(1-p))` remains for binaries, reported as
#' `inject_truth$<name>$binary_sim_factor`. Input ordered-discrete columns have
#' their latent correlations recovered so associations are preserved rather than
#' attenuated on the round-trip: binary-like columns (a numeric 0/1 column or a
#' two-level factor such as M/F) via tetrachoric / biserial, and ordered
#' multi-level columns (ordinal-mode integers, ordered factors) via polychoric /
#' polyserial. These recoveries assume a thresholded bivariate-normal latent.
#'
#' **Limitations.** Categorical phenotype injection is ordinal; for a nominal
#' grouping factor whose levels perturb different features independently, use
#' [add_group_effect()] (a post-processor). Unordered (>2-level) categorical
#' *input* columns are not polychoric-corrected (no inherent order), and
#' categorical input features are not supported as effect targets in this
#' version.
#'
#' @return An object of class `"synthetica_sim"` (a list) with elements:
#' * `data` -- the synthetic `data.frame` (n x p)
#' * `marginals` -- fitted marginals, one per column
#' * `corr` -- the projected latent correlation matrix
#' * `missingness` -- missingness fit (rates and, if applicable, latent copula)
#' * `inject_truth` -- per-phenotype bookkeeping (`NULL` when no phenotype is injected)
#' * `dropped_cols`, `dropped_near_constant`, `dropped_high_missing` --
#'    audit trail of guardrail actions
#' * `audit` -- one-line summary string of guardrail actions
#' * `call` -- the matched call
#'
#' @examples
#' set.seed(1)
#' input <- data.frame(
#'   x   = rnorm(200, 50, 10),
#'   y   = rlnorm(200, log(20), 0.5),
#'   sex = sample(c(0L, 1L), 200, replace = TRUE),
#'   grp = factor(sample(letters[1:3], 200, replace = TRUE))
#' )
#' input$x[sample(200, 30)] <- NA
#'
#' # Baseline simulation: synthetic copy
#' sim <- simulate_dataset(input, n = 500, seed = 42, verbose = FALSE)
#' str(sim$data)
#'
#' # Phenotype injection: a synthetic age column correlated with all features
#' sim2 <- simulate_dataset(
#'   input,
#'   n = 500, seed = 42, verbose = FALSE,
#'   inject = list(
#'     age = list(
#'       type = "quantitative", mean = 55, sd = 12,
#'       effect = list(kind = "mean_shift", size_sd = 0.2, on = "all")
#'     )
#'   )
#' )
#' sim2$inject_truth$age$mean_abs_real_cor
#'
#' @export
#' @importFrom MASS ginv mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom stats approx cor dnorm na.omit pnorm qnorm quantile rnorm runif sd setNames
#' @importFrom utils modifyList
simulate_dataset <- function(
  data,
  n                     = nrow(data),
  col_types             = NULL,
  marginal_resolution   = 128L,
  missingness           = c("preserve_pattern", "MCAR", "none"),
  inject                = NULL,
  copula_regularization = "nearPD",
  privacy               = list(min_cell_size      = 10L,
                               noise_on_corr      = 0,
                               drop_near_constant = TRUE,
                               max_missingness    = 0.90,
                               force              = FALSE),
  seed                  = NULL,
  verbose               = TRUE,
  stratify_by           = NULL,
  min_strata_n          = 30L,
  max_ordinal_levels    = 10L
) {
  # ---- input validation ----
  stopifnot("data must be a data.frame or matrix" = is.data.frame(data) || is.matrix(data))
  if (is.matrix(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  stopifnot("n must be a positive integer scalar" =
              is.numeric(n) && length(n) == 1L && n >= 1)
  n <- as.integer(n)
  stopifnot("max_ordinal_levels must be a positive integer scalar" =
              is.numeric(max_ordinal_levels) && length(max_ordinal_levels) == 1L &&
              max_ordinal_levels >= 1)
  max_ordinal_levels <- as.integer(max_ordinal_levels)
  missingness <- match.arg(missingness)

  # ---- stratify_by dispatch ----
  if (!is.null(stratify_by)) {
    stopifnot("stratify_by must be a single column name" =
                is.character(stratify_by) && length(stratify_by) == 1L)
    if (!stratify_by %in% names(data)) {
      stop("stratify_by column '", stratify_by, "' not found in data")
    }
    if (!is.null(inject)) {
      stop("stratify_by cannot be combined with inject in this version; ",
           "use them in separate calls")
    }
    stopifnot("min_strata_n must be a positive integer scalar" =
                is.numeric(min_strata_n) && length(min_strata_n) == 1L &&
                min_strata_n >= 1)
    min_strata_n <- as.integer(min_strata_n)
    return(.simulate_stratified(
      data                  = data,
      n                     = n,
      stratify_by           = stratify_by,
      min_strata_n          = min_strata_n,
      col_types             = col_types,
      marginal_resolution   = marginal_resolution,
      missingness           = missingness,
      copula_regularization = copula_regularization,
      privacy               = privacy,
      seed                  = seed,
      verbose               = verbose,
      max_ordinal_levels    = max_ordinal_levels
    ))
  }

  priv_def <- list(min_cell_size = 10L, noise_on_corr = 0,
                   drop_near_constant = TRUE, max_missingness = 0.90,
                   force = FALSE)
  privacy  <- modifyList(priv_def, privacy %||% list())

  call_obj         <- match.call()
  set_seed_locally <- function(s) if (!is.null(s)) set.seed(s)

  # ---- type detection & privacy guardrails ----
  types   <- .detect_col_types(data, override = col_types)
  guarded <- .apply_privacy_guardrails(
    df                 = data,
    types              = types,
    min_cell_size      = privacy$min_cell_size,
    drop_near_constant = privacy$drop_near_constant,
    max_missingness    = privacy$max_missingness,
    force              = privacy$force,
    verbose            = verbose
  )
  data  <- guarded$data
  types <- guarded$types
  if (ncol(data) < 2L) stop("after privacy guardrails fewer than 2 columns remain -- cannot fit a copula.")

  # ---- phenotype injection: attach injected phenotypes ----
  inject_truth   <- NULL
  inject_latents <- list()
  if (!is.null(inject)) {
    if (isTRUE(verbose)) message("attaching injected phenotype(s)")
    set_seed_locally(if (!is.null(seed)) seed + 100L else NULL)
    att            <- .attach_phenotypes(data, types, inject, verbose)
    data           <- att$data
    types          <- att$types
    inject_truth   <- att$inject_truth
    inject_latents <- att$latents
  }
  p <- ncol(data)

  # ---- fit marginals ----
  if (isTRUE(verbose)) message(sprintf(
    "fitting marginals on %d column(s) (resolution = %d)", p, marginal_resolution))
  marginals <- setNames(vector("list", p), names(data))
  for (j in seq_len(p)) {
    marginals[[j]] <- .fit_marginal(data[[j]], types[j],
                                    resolution         = marginal_resolution,
                                    min_cell_size      = privacy$min_cell_size,
                                    max_ordinal_levels = max_ordinal_levels)
  }

  # ---- PIT each column to latent z ----
  set_seed_locally(seed)
  Z <- matrix(NA_real_, nrow = nrow(data), ncol = p, dimnames = list(NULL, names(data)))
  for (j in seq_len(p)) Z[, j] <- .pit_to_latent(data[[j]], marginals[[j]])

  # Splice the true generating latent for injected phenotypes (thresholded-
  # normal): the copula then captures the designed correlation exactly instead
  # of re-imputing (and attenuating) it from the realized values.
  for (nm in names(inject_latents)) Z[, nm] <- inject_latents[[nm]]

  # ---- latent correlation ----
  if (isTRUE(verbose)) message("estimating latent correlation matrix")
  R <- .estimate_latent_corr(Z, noise_on_corr = privacy$noise_on_corr)
  # Recover faithful latent correlations for input ordered-discrete columns
  # (binary: tetrachoric/biserial; ordered multi-level: polychoric/polyserial),
  # undoing the thresholding attenuation. Re-projects to the nearest PD
  # correlation matrix internally if it changes anything.
  R <- .correct_discrete_corr(R, data, types, marginals, injected = names(inject_latents))

  # ---- fit missingness ----
  if (isTRUE(verbose)) message(sprintf("fitting missingness model (mode = %s)", missingness))
  miss_fit <- .fit_missingness(data, missingness)

  # ---- simulate ----
  if (isTRUE(verbose)) message(sprintf("sampling %d row(s) from copula", n))
  set_seed_locally(if (!is.null(seed)) seed + 1L else NULL)
  Z_sim <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = R)
  if (is.vector(Z_sim)) Z_sim <- matrix(Z_sim, nrow = 1L)

  sim <- as.data.frame(
    lapply(seq_len(p), function(j) .back_transform(Z_sim[, j], marginals[[j]])),
    stringsAsFactors = FALSE
  )
  names(sim) <- names(data)

  # ---- apply missingness ----
  set_seed_locally(if (!is.null(seed)) seed + 2L else NULL)
  mask <- .sample_missingness(miss_fit, n = n, p = p, col_names = names(data))
  for (j in seq_len(p)) sim[[j]][mask[, j]] <- NA

  structure(
    list(
      data                  = sim,
      marginals             = marginals,
      corr                  = R,
      missingness           = miss_fit,
      inject_truth          = inject_truth,
      dropped_cols          = guarded$dropped_cols,
      dropped_near_constant = guarded$dropped_near_constant,
      dropped_high_missing  = guarded$dropped_high_missing,
      audit                 = guarded$audit,
      call                  = call_obj
    ),
    class = "synthetica_sim"
  )
}
