# F3b: integer / whole-number column handling.

test_that("detection: whole vs fractional vs binary", {
  df <- data.frame(
    whole_dbl = c(1, 2, 3, 4, 5, 6),       # whole double, >2 distinct -> integer
    int_typed = 1:6,                        # integer-typed -> integer
    frac      = c(1.5, 2.2, 3.9, 4.1, 5, 6),# fractional -> quantitative
    twolev    = c(0, 1, 0, 1, 0, 1)         # 2 distinct -> binary
  )
  ty <- synthetica:::.detect_col_types(df)
  expect_equal(unname(ty["whole_dbl"]), "integer")
  expect_equal(unname(ty["int_typed"]), "integer")
  expect_equal(unname(ty["frac"]),      "quantitative")
  expect_equal(unname(ty["twolev"]),    "binary")
})

test_that("ordinal integer: exact level set preserved, integer output", {
  set.seed(1); n <- 800
  df <- data.frame(
    score = sample(1:5, n, replace = TRUE, prob = c(.1, .2, .4, .2, .1)),
    x     = rnorm(n)
  )
  sim <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE)
  expect_equal(sim$marginals$score$mode, "ordinal")
  expect_true(is.integer(sim$data$score))
  expect_true(all(sort(unique(sim$data$score)) %in% 1:5))
  # frequencies approximately preserved (absolute tolerance on proportions)
  pin  <- prop.table(table(factor(df$score,        levels = 1:5)))
  pout <- prop.table(table(factor(sim$data$score,  levels = 1:5)))
  expect_lt(max(abs(pin - pout)), 0.06)
})

test_that("grid integer: integer output within observed range", {
  set.seed(2); n <- 1000
  df <- data.frame(
    age = round(rnorm(n, 50, 15)),   # many distinct -> grid
    x   = rnorm(n)
  )
  sim <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE)
  expect_equal(sim$marginals$age$mode, "grid")
  expect_true(is.integer(sim$data$age))
  expect_gte(min(sim$data$age), min(df$age))
  expect_lte(max(sim$data$age), max(df$age))
})

test_that("min_cell_size forces a rare-level integer column to the grid path", {
  set.seed(3); n <- 600
  # values 7,8,9 each appear once -> below min_cell_size -> grid, not ordinal
  df <- data.frame(
    v = sample(c(rep(2L, n - 3L), 7L, 8L, 9L)),
    x = rnorm(n)
  )
  sim <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE)
  expect_equal(sim$marginals$v$mode, "grid")
  expect_true(is.integer(sim$data$v))
})

test_that("max_ordinal_levels boundary pushes wider supports to the grid", {
  set.seed(4); n <- 2000
  df <- data.frame(
    k = sample(1:8, n, replace = TRUE),   # 8 distinct
    x = rnorm(n)
  )
  ord  <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE,
                           max_ordinal_levels = 8)
  grid <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE,
                           max_ordinal_levels = 5)
  expect_equal(ord$marginals$k$mode,  "ordinal")
  expect_equal(grid$marginals$k$mode, "grid")
  expect_true(is.integer(ord$data$k))
  expect_true(is.integer(grid$data$k))
})

test_that("col_types override to quantitative keeps floats", {
  set.seed(5); n <- 800
  df <- data.frame(
    age = round(rnorm(n, 50, 15)),
    x   = rnorm(n)
  )
  sim <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE,
                          col_types = list(age = "quantitative"))
  expect_true(is.double(sim$data$age))
  expect_true(any(sim$data$age != round(sim$data$age)))   # genuinely fractional
})

test_that("whole-valued binary column returns integer", {
  set.seed(6); n <- 500
  df <- data.frame(
    flag = sample(c(0L, 1L), n, replace = TRUE),
    x    = rnorm(n)
  )
  sim <- simulate_dataset(df, n = n, seed = 42, verbose = FALSE)
  expect_true(is.integer(sim$data$flag))
  expect_true(all(sim$data$flag %in% c(0L, 1L)))
})

test_that("mtcars round-trip: count-like columns are integer, fractional stay double", {
  sim <- simulate_dataset(mtcars, n = nrow(mtcars), seed = 42, verbose = FALSE)
  for (col in c("cyl", "hp", "gear", "carb", "vs", "am")) {
    expect_true(is.integer(sim$data[[col]]), info = col)
  }
  for (col in c("mpg", "disp", "drat", "wt", "qsec")) {
    expect_true(is.double(sim$data[[col]]), info = col)
  }
  # NB: cyl/gear/carb each have a level below min_cell_size (e.g. cyl == 6
  # appears 7x), so they fall to grid mode -- integer output within range, not
  # exact original support. Confirm the in-range integer contract.
  expect_gte(min(sim$data$cyl), min(mtcars$cyl))
  expect_lte(max(sim$data$cyl), max(mtcars$cyl))
})
