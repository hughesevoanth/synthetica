#' Simulate Longitudinal / Repeated-Measures Data
#'
#' Generate privacy-preserving synthetic copies of longitudinal (panel) data --
#' repeated measures on the same subjects over time -- by identifying a subject
#' `id` column and a `time` column (e.g. `Chick` and `Time` in
#' `datasets::ChickWeight`).
#'
#' The data are reshaped to **wide** format (one row per subject, with a column
#' for each time-varying feature at each timepoint), passed through
#' [simulate_dataset()], and reshaped back to long. This turns within-subject
#' temporal autocorrelation into ordinary between-column correlation, so the
#' full copula machinery applies for free: per-timepoint marginals capture the
#' trajectory shape (e.g. growth), the missingness model reproduces dropout,
#' time-invariant covariates ride along as single columns, and integer /
#' discrete / polychoric handling all work.
#'
#' This handles ChickWeight-style panels, before/after RCTs (two timepoints),
#' crossover trials (columns per period), and modest biomarker panels. It does
#' **not** scale to high-dimensional longitudinal 'omics (a warning fires when
#' `features x timepoints` exceeds a few hundred columns), which needs a
#' different method.
#'
#' @param data A long-format `data.frame` with one row per (subject, time).
#' @param id Character, the name of the subject identifier column.
#' @param time Character, the name of the time column. Its unique observed
#'   values define the time grid; genuinely irregular/continuous times should be
#'   binned to a common grid before calling.
#' @param n Number of synthetic **subjects** to generate. Default `NULL` uses
#'   the number of input subjects.
#' @param static Optional character vector of time-invariant covariate columns
#'   (constant within each subject, e.g. `Diet`). `NULL` (default) auto-detects
#'   them. Columns declared static that actually vary within a subject emit a
#'   warning and use the first value per subject.
#' @param id_prefix Character prefix for the synthetic subject ids. Default
#'   `"subj"` (ids become `subj1`, `subj2`, ...).
#' @param verbose Logical; passed to [simulate_dataset()].
#' @param ... Further arguments forwarded to [simulate_dataset()] (`seed`,
#'   `privacy`, `missingness`, `marginal_resolution`, `max_ordinal_levels`, ...).
#'
#' @return An object of class `"synthetica_long"`: a list with
#' * `data` -- the synthetic long-format `data.frame` (same columns as input)
#' * `wide` -- the underlying `synthetica_sim` fitted on the wide frame
#' * `id`, `time`, `times`, `static`, `features`, `n_subjects` -- metadata
#'
#' @examples
#' data(ChickWeight, package = "datasets")
#' sim <- simulate_longitudinal(ChickWeight, id = "Chick", time = "Time",
#'                              seed = 1, verbose = FALSE)
#' head(sim$data)
#' # synthetic chicks have growing weight trajectories and a Diet covariate
#'
#' @seealso [simulate_dataset()]
#' @export
simulate_longitudinal <- function(data, id, time, n = NULL, static = NULL,
                                  id_prefix = "subj", verbose = TRUE, ...) {
  stopifnot("data must be a data.frame" = is.data.frame(data))
  stopifnot("id must be a single column name present in data" =
              is.character(id) && length(id) == 1L && id %in% names(data))
  stopifnot("time must be a single column name present in data" =
              is.character(time) && length(time) == 1L && time %in% names(data))
  stopifnot("id and time must differ" = id != time)

  ids   <- unique(data[[id]])
  times <- sort(unique(data[[time]]))
  n_subj <- length(ids)
  if (is.null(n)) n <- n_subj
  stopifnot("n must be a positive integer scalar" =
              is.numeric(n) && length(n) == 1L && n >= 1)
  n <- as.integer(n)
  if (length(times) < 2L) stop("need at least 2 distinct time points")

  # ---- column roles ----
  other <- setdiff(names(data), c(id, time))
  if (!length(other)) stop("no feature columns besides id and time")
  const_within_id <- vapply(other, function(col) {
    all(tapply(data[[col]], data[[id]],
               function(v) length(unique(v[!is.na(v)])) <= 1L))
  }, logical(1L))
  detected_static <- other[const_within_id]

  if (is.null(static)) {
    static <- detected_static
  } else {
    miss <- setdiff(static, other)
    if (length(miss)) stop("static column(s) not in data: ", paste(miss, collapse = ", "))
    varying <- setdiff(static, detected_static)
    if (length(varying)) {
      warning("declared static but varies within id (using first value per subject): ",
              paste(varying, collapse = ", "))
    }
  }
  features <- setdiff(other, static)
  if (!length(features)) stop("no time-varying feature columns (all are static?)")

  n_wide <- length(features) * length(times)
  if (n_wide > 500L) {
    warning(sprintf(
      "%d feature x time columns: wide-format is slow/infeasible above a few hundred; high-dimensional longitudinal 'omics is out of scope for this function.",
      n_wide))
  }

  # ---- pivot to wide (one row per subject) ----
  cmap <- expand.grid(feature = features, time = times, stringsAsFactors = FALSE)
  cmap$wide_col <- make.unique(make.names(paste0(cmap$feature, "_T", cmap$time)))

  first_idx <- match(ids, data[[id]])
  wide <- data.frame(row.names = seq_len(n_subj))
  for (r in seq_len(nrow(cmap))) {
    sub <- data[data[[time]] == cmap$time[r], , drop = FALSE]
    wide[[cmap$wide_col[r]]] <- sub[[cmap$feature[r]]][match(ids, sub[[id]])]
  }
  for (s in static) wide[[s]] <- data[[s]][first_idx]

  # ---- simulate on the wide frame ----
  if (isTRUE(verbose)) {
    message(sprintf("longitudinal: %d subjects x %d time(s); %d time-varying feature(s), %d static; %d wide column(s)",
                    n_subj, length(times), length(features), length(static), ncol(wide)))
  }
  wide_sim <- simulate_dataset(wide, n = n, verbose = verbose, ...)
  syn_wide <- wide_sim$data

  # ---- pivot back to long ----
  syn_ids <- factor(paste0(id_prefix, seq_len(n)),
                    levels = paste0(id_prefix, seq_len(n)))
  blocks <- lapply(times, function(t) {
    blk <- stats::setNames(
      data.frame(syn_ids, t, stringsAsFactors = FALSE, check.names = FALSE),
      c(id, time))
    for (f in features) {
      wc <- cmap$wide_col[cmap$feature == f & cmap$time == t]
      blk[[f]] <- if (length(wc) && wc %in% names(syn_wide)) syn_wide[[wc]] else NA
    }
    for (s in static) blk[[s]] <- if (s %in% names(syn_wide)) syn_wide[[s]] else NA
    blk
  })
  long <- do.call(rbind, blocks)

  # drop (subject, time) rows where no time-varying feature was observed (dropout)
  keep <- rowSums(!is.na(long[, features, drop = FALSE])) > 0
  long <- long[keep, , drop = FALSE]
  long <- long[, intersect(names(data), names(long)), drop = FALSE]
  # order by subject then time for readability
  long <- long[order(long[[id]], long[[time]]), , drop = FALSE]
  rownames(long) <- NULL

  structure(
    list(data = long, wide = wide_sim,
         id = id, time = time, times = times,
         static = static, features = features, n_subjects = n),
    class = "synthetica_long")
}
