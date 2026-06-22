# F3a: thresholded-normal binary fidelity -- helper unit tests and round-trip
# preservation. Internal (dotted) helpers are reachable under devtools::test().

test_that(".tetrachoric_digby recovers a known latent correlation", {
  set.seed(1)
  n <- 5000; rho <- 0.5
  z1 <- rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  x  <- as.integer(z1 > 0)                  # prevalence ~0.5
  y  <- as.integer(z2 > qnorm(0.7))         # prevalence ~0.3
  est <- .tetrachoric_digby(x, y)
  expect_lt(abs(est - rho), 0.10)
})

test_that(".polyserial_latent recovers a known latent correlation (binary case = biserial)", {
  set.seed(2)
  n <- 5000; rho <- 0.5
  z1 <- rnorm(n)
  cont <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  b    <- as.integer(z1 > qnorm(0.6))       # prevalence ~0.4
  est  <- .polyserial_latent(b, cont)
  expect_lt(abs(est - rho), 0.08)
})

test_that("binary helpers clip / stay finite on degenerate inputs", {
  # empty cell (continuity correction keeps it finite)
  x <- c(rep(1L, 60), rep(0L, 40)); y <- c(rep(1L, 60), rep(0L, 40))
  expect_true(is.finite(.tetrachoric_digby(x, y)))
  expect_lte(abs(.tetrachoric_digby(x, y)), 0.999)
  # extreme prevalence polyserial stays finite and clipped
  b <- c(rep(1L, 2), rep(0L, 98)); cc <- rnorm(100)
  est <- .polyserial_latent(b, cc)
  expect_true(is.finite(est))
  expect_lte(abs(est), 0.999)
})

test_that("Baseline:input binary-continuous correlation is preserved, not attenuated", {
  set.seed(3)
  n <- 2000; rho <- 0.6
  z  <- rnorm(n)
  df <- data.frame(
    cont = rho * z + sqrt(1 - rho^2) * rnorm(n),
    bin  = as.integer(z > 0)
  )
  obs <- cor(df$bin, df$cont)                        # input point-biserial
  sim <- simulate_dataset(df, n = 4000, seed = 42, verbose = FALSE)
  syn <- cor(as.numeric(sim$data$bin), sim$data$cont)
  expect_lt(abs(syn - obs), 0.06)
})

test_that("Baseline:input binary-binary correlation is preserved", {
  set.seed(4)
  n <- 3000; rho <- 0.5
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  df <- data.frame(
    a = as.integer(z1 > 0),
    b = as.integer(z2 > qnorm(0.6)),
    c = rnorm(n)
  )
  obs <- cor(df$a, df$b)
  sim <- simulate_dataset(df, n = 5000, seed = 42, verbose = FALSE)
  syn <- cor(as.numeric(sim$data$a), as.numeric(sim$data$b))
  expect_lt(abs(syn - obs), 0.06)
})

test_that("Baseline:two-level factor (M/F) is treated as binary-like and preserved", {
  set.seed(6)
  n <- 3000; rho <- 0.55
  z  <- rnorm(n)
  df <- data.frame(
    cont = rho * z + sqrt(1 - rho^2) * rnorm(n),
    sex  = factor(ifelse(z > 0, "M", "F"))
  )
  obs <- cor(as.integer(df$sex), df$cont)
  sim <- simulate_dataset(df, n = 4000, seed = 42, verbose = FALSE)
  expect_s3_class(sim$data$sex, "factor")
  syn <- cor(as.integer(sim$data$sex), sim$data$cont)
  expect_lt(abs(syn - obs), 0.06)
})

test_that("Phenotype injection: injected binary realizes correlation near the predicted single factor", {
  set.seed(5)
  n  <- 2000
  df <- data.frame(x = rnorm(n), y = rnorm(n))
  sim <- simulate_dataset(
    df, n = 4000, seed = 42, verbose = FALSE,
    inject = list(
      g = list(type = "binary", prob = 0.5,
               effect = list(kind = "r_squared", min_r2 = 0.3, on = c("x", "y")))
    )
  )
  pred <- sim$inject_truth$g$expected_sim_cor       # latent_target_cor * phi/sqrt(pq)
  syn  <- mean(abs(c(cor(as.numeric(sim$data$g), sim$data$x),
                     cor(as.numeric(sim$data$g), sim$data$y))))
  expect_lt(abs(syn - pred), 0.10)
})
