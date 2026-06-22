test_that("stratify_by returns a synthetica_sim with per-level marginals and corr", {
  sim <- simulate_dataset(iris, n = 300, seed = 42, verbose = FALSE,
                          stratify_by = "Species")
  expect_s3_class(sim, "synthetica_sim")
  expect_equal(nrow(sim$data), 300L)
  expect_true("Species" %in% names(sim$data))
  expect_setequal(levels(sim$data$Species), levels(iris$Species))
  expect_named(sim$marginals, levels(iris$Species))
  expect_named(sim$corr,      levels(iris$Species))
  # Per-stratum corr matrices: 4 numeric columns -> 4x4 each
  for (lev in levels(iris$Species)) {
    expect_equal(dim(sim$corr[[lev]]), c(4L, 4L))
  }
})

test_that("stratify_by preserves within-stratum medians", {
  sim <- simulate_dataset(iris, n = 3000, seed = 42, verbose = FALSE,
                          stratify_by = "Species")
  for (sp in levels(iris$Species)) {
    for (col in c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")) {
      in_med  <- median(iris[[col]][iris$Species == sp])
      sim_med <- median(sim$data[[col]][sim$data$Species == sp])
      expect_lt(abs(in_med - sim_med), 0.3,
                label = sprintf("%s within %s", col, sp))
    }
  }
})

test_that("stratify_by preserves per-level proportions", {
  sim <- simulate_dataset(iris, n = 3000, seed = 42, verbose = FALSE,
                          stratify_by = "Species")
  in_prop  <- prop.table(table(iris$Species))
  sim_prop <- prop.table(table(sim$data$Species))
  expect_lt(max(abs(in_prop - sim_prop)), 0.02)
})

test_that("stratify_by + inject errors in v1", {
  expect_error(
    simulate_dataset(
      iris, n = 100, verbose = FALSE,
      stratify_by = "Species",
      inject = list(z = list(type = "quantitative", mean = 0, sd = 1,
                             effect = list(kind = "mean_shift", size_sd = 0.1,
                                           on = "all")))
    ),
    "stratify_by cannot be combined with inject"
  )
})

test_that("stratify_by: small strata are skipped with a warning", {
  set.seed(1)
  input <- rbind(
    data.frame(x = rnorm(80, 0,  1), y = rnorm(80, 0,  1), g = "A"),
    data.frame(x = rnorm(80, 5,  1), y = rnorm(80, 5,  1), g = "B"),
    data.frame(x = rnorm(5,  10, 1), y = rnorm(5,  10, 1), g = "C")
  )
  input$g <- factor(input$g)
  expect_warning(
    sim <- simulate_dataset(input, n = 200, seed = 1, verbose = FALSE,
                            stratify_by = "g", min_strata_n = 30L),
    "strata with < 30 rows are skipped"
  )
  expect_false("C" %in% as.character(sim$data$g))
  expect_true(all(c("A", "B") %in% as.character(sim$data$g)))
  expect_equal(sim$small_levels_skipped, "C")
})

test_that("stratify_by: errors on missing column name", {
  expect_error(
    simulate_dataset(iris, n = 100, verbose = FALSE,
                     stratify_by = "NotAColumn"),
    "not found in data"
  )
})

test_that("stratify_by: NA values in stratify column are dropped with a warning", {
  input <- iris
  input$Species[c(10L, 60L, 110L)] <- NA
  expect_warning(
    sim <- simulate_dataset(input, n = 200, seed = 1, verbose = FALSE,
                            stratify_by = "Species"),
    "have NA in stratify column"
  )
  expect_false(anyNA(sim$data$Species))
})

test_that("stratify_by: output rows stay grouped by stratum", {
  sim <- simulate_dataset(iris, n = 300, seed = 42, verbose = FALSE,
                          stratify_by = "Species")
  rl <- rle(as.character(sim$data$Species))
  # Each level should appear as a single contiguous run
  expect_equal(length(rl$lengths), length(unique(rl$values)))
})

test_that("stratify_by: column order matches input", {
  sim <- simulate_dataset(iris, n = 100, seed = 1, verbose = FALSE,
                          stratify_by = "Species")
  expect_equal(names(sim$data), names(iris))
})
