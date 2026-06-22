# Polychoric / polyserial correction for ordered multi-level discrete columns.

cut_ord <- function(z, probs) {
  # cut a latent into ordinal codes 1..k with the given marginal proportions
  tau <- qnorm(cumsum(probs))[-length(probs)]
  findInterval(z, tau) + 1L
}

test_that(".bvn_cdf matches the closed form Phi2(0,0;rho) and the rho limits", {
  for (rho in c(-0.8, -0.3, 0, 0.3, 0.8)) {
    expect_lt(abs(.bvn_cdf(0, 0, rho) - (0.25 + asin(rho) / (2 * pi))), 1e-6)
  }
  # rho = 0 -> independence
  expect_lt(abs(.bvn_cdf(0.4, -0.7, 0) - pnorm(0.4) * pnorm(-0.7)), 1e-10)
  # rho -> 1 -> comonotone: Phi2(h,k;1) = min(Phi(h), Phi(k))
  expect_lt(abs(.bvn_cdf(0.4, -0.7, 0.999) - min(pnorm(0.4), pnorm(-0.7))), 0.01)
})

test_that(".polyserial_latent recovers a known latent correlation (k = 4)", {
  set.seed(1)
  n <- 6000; rho <- 0.55
  z   <- rnorm(n)
  cont <- rho * z + sqrt(1 - rho^2) * rnorm(n)
  ord  <- cut_ord(z, c(0.25, 0.25, 0.25, 0.25))
  est  <- .polyserial_latent(ord, cont)
  expect_lt(abs(est - rho), 0.06)
})

test_that(".polychoric_latent recovers a known latent correlation (k = 4 x m = 3)", {
  set.seed(2)
  n <- 8000; rho <- 0.5
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  a  <- cut_ord(z1, c(0.2, 0.3, 0.3, 0.2))
  b  <- cut_ord(z2, c(0.4, 0.3, 0.3))
  est <- .polychoric_latent(a, b)
  expect_lt(abs(est - rho), 0.06)
})

test_that(".polychoric_latent generalises tetrachoric at k = m = 2", {
  set.seed(3)
  n <- 6000; rho <- 0.45
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  a <- as.integer(z1 > 0); b <- as.integer(z2 > qnorm(0.6))
  expect_lt(abs(.polychoric_latent(a, b) - rho), 0.07)
})

test_that("Baseline:ordinal-integer <-> continuous correlation is preserved", {
  set.seed(4)
  n <- 3000; rho <- 0.6
  z  <- rnorm(n)
  df <- data.frame(
    score = cut_ord(z, rep(0.2, 5)),                  # 5-level ordinal integer
    y     = rho * z + sqrt(1 - rho^2) * rnorm(n)
  )
  # score has 5 well-populated levels -> ordinal mode
  obs <- cor(df$score, df$y, method = "spearman")
  sim <- simulate_dataset(df, n = 6000, seed = 42, verbose = FALSE)
  expect_equal(sim$marginals$score$mode, "ordinal")
  syn <- cor(as.numeric(sim$data$score), sim$data$y, method = "spearman")
  expect_lt(abs(syn - obs), 0.06)
})

test_that("Baseline:ordinal <-> ordinal correlation is preserved", {
  set.seed(5)
  n <- 4000; rho <- 0.55
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  df <- data.frame(
    s1 = cut_ord(z1, rep(0.25, 4)),
    s2 = cut_ord(z2, rep(0.2, 5)),
    x  = rnorm(n)
  )
  obs <- cor(df$s1, df$s2, method = "spearman")
  sim <- simulate_dataset(df, n = 8000, seed = 42, verbose = FALSE)
  expect_equal(sim$marginals$s1$mode, "ordinal")
  expect_equal(sim$marginals$s2$mode, "ordinal")
  syn <- cor(as.numeric(sim$data$s1), as.numeric(sim$data$s2), method = "spearman")
  expect_lt(abs(syn - obs), 0.07)
})
