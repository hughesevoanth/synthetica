test_that("Phenotype injection:injected quantitative phenotype appears in output", {
  input <- make_test_input(n = 500)
  sim <- simulate_dataset(
    input, n = 500, seed = 42, verbose = FALSE,
    inject = list(
      age = list(
        type = "quantitative", mean = 55, sd = 12,
        effect = list(kind = "mean_shift", size_sd = 0.20, on = "all")
      )
    )
  )
  expect_true("age" %in% names(sim$data))
  expect_equal(mean(sim$data$age, na.rm = TRUE), 55, tolerance = 2)
  expect_equal(sd(sim$data$age,   na.rm = TRUE), 12, tolerance = 2)
  expect_true("age" %in% names(sim$inject_truth))
  expect_equal(sim$inject_truth$age$binary_sim_factor, 1)
})

test_that("Phenotype injection:r_squared spec yields realized |cor| near sqrt(min_r2)", {
  set.seed(11)
  n <- 1000
  z <- rnorm(n)
  input <- data.frame(
    a = 0.3 * z + sqrt(1 - 0.09) * rnorm(n),
    b = 0.3 * z + sqrt(1 - 0.09) * rnorm(n),
    c = 0.3 * z + sqrt(1 - 0.09) * rnorm(n)
  )
  sim <- simulate_dataset(
    input, n = 1000, seed = 42, verbose = FALSE,
    inject = list(
      bmi = list(
        type = "quantitative", mean = 25, sd = 5,
        effect = list(kind = "r_squared", min_r2 = 0.10, on = "all")
      )
    )
  )
  # target |corr| = sqrt(0.10) ~ 0.316; expect within 30% (sample noise + ginv)
  expect_gt(sim$inject_truth$bmi$mean_abs_real_cor, 0.20)
})

test_that("Phenotype injection:binary phenotype has correct prevalence and thresholding factor", {
  input <- make_test_input(n = 1000)
  sim <- simulate_dataset(
    input, n = 1000, seed = 42, verbose = FALSE,
    inject = list(
      group = list(
        type = "binary", prob = 0.40,
        effect = list(kind = "mean_shift", size_sd = 0.10, on = "all")
      )
    )
  )
  # Absolute tolerance: binomial sampling SD for n=1000, p=0.4 is ~0.015,
  # so allow up to ~3 SDs.
  expect_lt(abs(mean(sim$data$group == 1, na.rm = TRUE) - 0.40), 0.05)
  # Single (unavoidable) thresholding factor phi(c) / sqrt(p*(1-p)) for p = 0.4,
  # now that the generating latent is spliced (no extra squared attenuation).
  p <- 0.40
  phi <- dnorm(qnorm(1 - p))
  expected <- phi / sqrt(p * (1 - p))
  expect_equal(sim$inject_truth$group$binary_sim_factor, expected, tolerance = 1e-6)
})

test_that("Phenotype injection:phenotype name colliding with feature errors", {
  input <- make_test_input()
  expect_error(
    simulate_dataset(
      input, n = 100, verbose = FALSE,
      inject = list(
        x = list(
          type = "quantitative", mean = 0, sd = 1,
          effect = list(kind = "mean_shift", size_sd = 0.1, on = "all")
        )
      )
    ),
    "collide"
  )
})

test_that("Phenotype injection:infeasible effect spec emits clipping warning", {
  input <- make_test_input(n = 500)
  expect_warning(
    simulate_dataset(
      input, n = 100, seed = 42, verbose = FALSE,
      inject = list(
        big = list(
          type = "quantitative", mean = 0, sd = 1,
          effect = list(kind = "mean_shift", size_sd = 0.95, on = "all")
        )
      )
    ),
    "infeasible"
  )
})

test_that("Phenotype injection:bad spec types are rejected", {
  input <- make_test_input()
  expect_error(
    simulate_dataset(
      input, n = 50, verbose = FALSE,
      inject = list(z = list(type = "ordinal", mean = 0, sd = 1,
                             effect = list(kind = "mean_shift", size_sd = 0.1, on = "all")))
    ),
    "type"
  )
  expect_error(
    simulate_dataset(
      input, n = 50, verbose = FALSE,
      inject = list(z = list(type = "binary", prob = 1.5,
                             effect = list(kind = "mean_shift", size_sd = 0.1, on = "all")))
    ),
    "prob"
  )
})

# ---- Model A: ordinal categorical injection ------------------------------

test_that("Phenotype injection:injected categorical phenotype appears as a factor with right levels/proportions", {
  input <- make_test_input(n = 1000)
  sim <- simulate_dataset(
    input, n = 1000, seed = 42, verbose = FALSE,
    inject = list(
      soil = list(
        type = "categorical", levels = c("clay", "sand", "silt"),
        probs = c(0.5, 0.3, 0.2),
        effect = list(kind = "r_squared", min_r2 = 0.3, on = c("x", "y"))
      )
    )
  )
  expect_true("soil" %in% names(sim$data))
  expect_s3_class(sim$data$soil, "factor")
  expect_equal(levels(sim$data$soil), c("clay", "sand", "silt"))
  props <- as.numeric(prop.table(table(sim$data$soil)))
  expect_lt(max(abs(props - c(0.5, 0.3, 0.2))), 0.05)
  expect_equal(sim$inject_truth$soil$binary_sim_factor, 1)
})

test_that("Phenotype injection:categorical phenotype is ordinally correlated with its targets", {
  input <- make_test_input(n = 1000)
  sim <- simulate_dataset(
    input, n = 1000, seed = 42, verbose = FALSE,
    inject = list(
      soil = list(
        type = "categorical", levels = c("clay", "sand", "silt"),
        probs = c(0.4, 0.3, 0.3),
        effect = list(kind = "r_squared", min_r2 = 0.4, on = c("x", "y"))
      )
    )
  )
  # realized association recorded per target, in the designed (positive) direction
  rc <- sim$inject_truth$soil$realized_cor
  expect_named(rc, c("x", "y"))
  expect_gt(min(rc), 0.2)
  # group medians should increase along the level order (clay < sand < silt)
  meds <- tapply(sim$data$x, sim$data$soil, median)
  expect_true(meds[["clay"]] < meds[["sand"]] && meds[["sand"]] < meds[["silt"]])
})

test_that("Phenotype injection:stronger min_r2 yields stronger realized categorical association", {
  input <- make_test_input(n = 1000)
  mk <- function(r2) simulate_dataset(
    input, n = 1000, seed = 42, verbose = FALSE,
    inject = list(soil = list(
      type = "categorical", levels = c("a", "b", "c"), probs = c(1/3, 1/3, 1/3),
      effect = list(kind = "r_squared", min_r2 = r2, on = c("x", "y"))
    ))
  )$inject_truth$soil$mean_abs_real_cor
  expect_gt(mk(0.3), mk(0.05))
})

test_that("Phenotype injection:bad categorical specs are rejected", {
  input <- make_test_input()
  base_eff <- list(kind = "r_squared", min_r2 = 0.3, on = "x")

  # probs do not sum to 1
  expect_error(
    simulate_dataset(input, n = 50, verbose = FALSE, inject = list(
      s = list(type = "categorical", levels = c("a", "b"), probs = c(0.5, 0.4),
               effect = base_eff))),
    "sum to 1"
  )
  # levels / probs length mismatch
  expect_error(
    simulate_dataset(input, n = 50, verbose = FALSE, inject = list(
      s = list(type = "categorical", levels = c("a", "b", "c"), probs = c(0.5, 0.5),
               effect = base_eff))),
    "same length"
  )
  # duplicate levels
  expect_error(
    simulate_dataset(input, n = 50, verbose = FALSE, inject = list(
      s = list(type = "categorical", levels = c("a", "a"), probs = c(0.5, 0.5),
               effect = base_eff))),
    "unique"
  )
  # mean_shift not allowed for categorical
  expect_error(
    simulate_dataset(input, n = 50, verbose = FALSE, inject = list(
      s = list(type = "categorical", levels = c("a", "b"), probs = c(0.5, 0.5),
               effect = list(kind = "mean_shift", size_sd = 0.3, on = "x")))),
    "r_squared"
  )
})
