test_that("Privacy: drops columns with missingness above max_missingness", {
  input <- make_test_input(n = 100)
  input$mostly_na <- NA_real_
  input$mostly_na[1:5] <- rnorm(5)  # 95% missing > default 90%
  sim <- simulate_dataset(input, n = 50, seed = 42, verbose = FALSE)
  expect_false("mostly_na" %in% names(sim$data))
  expect_true("mostly_na" %in% sim$dropped_high_missing)
})

test_that("Privacy: drops near-constant columns", {
  input <- make_test_input(n = 100)
  input$const <- 5
  sim <- simulate_dataset(input, n = 50, seed = 42, verbose = FALSE)
  expect_false("const" %in% names(sim$data))
  expect_true("const" %in% sim$dropped_near_constant)
})

test_that("Privacy: collapses rare categorical levels into 'other'", {
  set.seed(7)
  input <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    g = factor(c(rep("a", 80), rep("b", 18), rep("c", 2)))
  )
  sim <- simulate_dataset(
    input, n = 50, seed = 42, verbose = FALSE,
    privacy = list(min_cell_size = 5)
  )
  expect_true("other" %in% levels(sim$data$g))
  expect_false("c"     %in% levels(sim$data$g))
})

test_that("Privacy: min_cell_size < 5 errors without force", {
  input <- make_test_input()
  expect_error(
    simulate_dataset(
      input, n = 50, verbose = FALSE,
      privacy = list(min_cell_size = 2)
    ),
    "force"
  )
})

test_that("Privacy: min_cell_size < 5 with force = TRUE warns and runs", {
  input <- make_test_input()
  expect_warning(
    sim <- simulate_dataset(
      input, n = 50, verbose = FALSE,
      privacy = list(min_cell_size = 2, force = TRUE)
    ),
    "re-identification"
  )
  expect_s3_class(sim, "synthetica_sim")
})

test_that("Privacy: marginals store quantile grids, not raw values", {
  input <- make_test_input(n = 500)
  sim   <- simulate_dataset(input, n = 100, seed = 42, verbose = FALSE,
                            marginal_resolution = 64L)
  m <- sim$marginals$x
  expect_equal(m$type, "quantitative")
  expect_equal(length(m$quant_grid), 65L)
  expect_equal(length(m$prob_grid),  65L)
  # Stored quantile grid should NOT equal the raw input vector
  expect_false(length(m$quant_grid) == length(input$x))
})
