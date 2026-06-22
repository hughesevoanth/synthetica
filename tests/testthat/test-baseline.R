test_that("simulate_dataset returns expected shape", {
  input <- make_test_input()
  sim   <- simulate_dataset(input, n = 100, seed = 42, verbose = FALSE)
  expect_s3_class(sim, "synthetica_sim")
  expect_equal(nrow(sim$data), 100L)
  expect_setequal(names(sim$data), names(input))
  expect_true(is.list(sim$marginals))
  expect_true(is.matrix(sim$corr))
})

test_that("Baseline simulation preservesmarginal medians within tolerance", {
  input <- make_test_input(n = 1000)
  sim   <- simulate_dataset(input, n = 5000, seed = 42, verbose = FALSE)
  for (nm in c("x", "y")) {
    expect_equal(
      median(sim$data[[nm]], na.rm = TRUE),
      median(input[[nm]],    na.rm = TRUE),
      tolerance = 0.10, info = nm
    )
  }
  expect_equal(
    mean(sim$data$sex == 1, na.rm = TRUE),
    mean(input$sex == 1, na.rm = TRUE),
    tolerance = 0.05
  )
})

test_that("Baseline simulation preservesSpearman correlation between correlated features", {
  set.seed(2)
  n <- 1000
  x <- rnorm(n)
  y <- 0.6 * x + sqrt(1 - 0.36) * rnorm(n)
  input <- data.frame(x = x, y = y)
  sim   <- simulate_dataset(input, n = 5000, seed = 42, verbose = FALSE)
  in_rho  <- cor(input$x,    input$y,    method = "spearman")
  sim_rho <- cor(sim$data$x, sim$data$y, method = "spearman")
  # Absolute tolerance: nearPD projection + rank-transform discretization can
  # shrink the realized cor by ~0.04 even at large n. 0.10 is the acceptance band.
  expect_lt(abs(sim_rho - in_rho), 0.10)
})

test_that("Baseline simulation propagates per-column missingness rate", {
  input <- make_test_input(n = 500)
  input$x[sample(500, 100)] <- NA
  sim <- simulate_dataset(input, n = 1000, seed = 42, verbose = FALSE)
  expect_lt(abs(mean(is.na(sim$data$x)) - mean(is.na(input$x))), 0.05)
})

test_that("Baseline:no exact-row match between input and sim", {
  input <- make_numeric_fixture(n = 200)
  sim   <- simulate_dataset(input, n = 200, seed = 42, verbose = FALSE)
  hash <- function(M) apply(format(M, digits = 12), 1L, paste, collapse = "|")
  in_hash  <- hash(as.matrix(input))
  sim_hash <- hash(as.matrix(sim$data))
  expect_equal(sum(sim_hash %in% in_hash), 0L)
})

test_that("Baseline:seed produces reproducible output", {
  input <- make_test_input()
  s1 <- simulate_dataset(input, n = 100, seed = 7, verbose = FALSE)
  s2 <- simulate_dataset(input, n = 100, seed = 7, verbose = FALSE)
  expect_equal(s1$data, s2$data)
})

test_that("Baseline:matrix input is accepted", {
  set.seed(3)
  M <- matrix(rnorm(400), nrow = 100, ncol = 4,
              dimnames = list(NULL, paste0("v", 1:4)))
  sim <- simulate_dataset(M, n = 50, seed = 1, verbose = FALSE)
  expect_s3_class(sim, "synthetica_sim")
  expect_equal(nrow(sim$data), 50L)
})
