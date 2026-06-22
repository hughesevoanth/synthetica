# simulate_longitudinal(): wide-format repeated-measures simulation.

# correlation of a feature between two timepoints, aligned by subject
pair_cor <- function(df, id, time, val, t0, t1) {
  a <- df[df[[time]] == t0, c(id, val)]
  b <- df[df[[time]] == t1, c(id, val)]
  m <- merge(a, b, by = id)
  cor(m[[2]], m[[3]], use = "complete.obs")
}

test_that("ChickWeight: structure, growth, trajectory correlation, static covariate", {
  data(ChickWeight, package = "datasets")
  sim <- simulate_longitudinal(ChickWeight, id = "Chick", time = "Time",
                               seed = 1, verbose = FALSE)
  expect_s3_class(sim, "synthetica_long")
  expect_equal(names(sim$data), names(ChickWeight))
  expect_equal(sim$static, "Diet")
  expect_equal(sim$features, "weight")
  expect_equal(length(unique(sim$data$Chick)), 50L)

  # growth: synthetic mean weight increases over time, like the real data
  mt <- tapply(sim$data$weight, sim$data$Time, mean)
  expect_gt(mt[["21"]], mt[["0"]])
  # per-time means track the real means (absolute tolerance in grams)
  rt <- tapply(ChickWeight$weight, ChickWeight$Time, mean)
  expect_lt(max(abs(mt[names(rt)] - rt)), 15)

  # within-subject trajectory correlation preserved
  rc <- pair_cor(ChickWeight, "Chick", "Time", "weight", 0, 20)
  sc <- pair_cor(sim$data,    "Chick", "Time", "weight", 0, 20)
  expect_lt(abs(sc - rc), 0.15)

  # Diet carried through as a factor with the same levels
  expect_s3_class(sim$data$Diet, "factor")
  expect_setequal(levels(sim$data$Diet), levels(ChickWeight$Diet))
})

test_that("two-timepoint RCT: pre/post correlation is preserved", {
  set.seed(1); ns <- 400
  base <- rnorm(ns, 10, 3)
  post <- 0.7 * base + sqrt(1 - 0.49) * rnorm(ns, 0, 3) + 5   # treatment effect + corr
  rct <- data.frame(
    subj  = rep(seq_len(ns), each = 2),
    visit = rep(c(0L, 1L), ns),
    y     = as.vector(rbind(base, post))
  )
  sim <- simulate_longitudinal(rct, id = "subj", time = "visit",
                               seed = 1, verbose = FALSE)
  rc <- pair_cor(rct,      "subj", "visit", "y", 0, 1)
  sc <- pair_cor(sim$data, "subj", "visit", "y", 0, 1)
  expect_lt(abs(sc - rc), 0.08)
  # post mean > pre mean (the treatment shift) preserved
  expect_gt(mean(sim$data$y[sim$data$visit == 1]),
            mean(sim$data$y[sim$data$visit == 0]))
})

test_that("static auto-detection and n control", {
  set.seed(2); ns <- 60
  df <- data.frame(
    id  = rep(seq_len(ns), each = 3),
    t   = rep(1:3, ns),
    x   = rnorm(ns * 3),
    grp = rep(sample(c("a", "b"), ns, TRUE), each = 3)   # constant within id
  )
  sim <- simulate_longitudinal(df, id = "id", time = "t", n = 20,
                               seed = 1, verbose = FALSE)
  expect_equal(sim$static, "grp")
  expect_equal(sim$features, "x")
  expect_equal(length(unique(sim$data$id)), 20L)
})

test_that("declaring a time-varying column as static warns", {
  set.seed(3); ns <- 50
  df <- data.frame(
    id = rep(seq_len(ns), each = 3),
    t  = rep(1:3, ns),
    x  = rnorm(ns * 3),                                   # varies within id
    z  = rnorm(ns * 3)
  )
  expect_warning(
    simulate_longitudinal(df, id = "id", time = "t", static = "x",
                          seed = 1, verbose = FALSE),
    "varies within id"
  )
})

test_that("errors on missing id/time and too-few timepoints", {
  df <- data.frame(id = 1:10, t = 1, y = rnorm(10))
  expect_error(simulate_longitudinal(df, id = "nope", time = "t"), "id")
  expect_error(simulate_longitudinal(df, id = "id", time = "t"), "2 distinct time")
})
