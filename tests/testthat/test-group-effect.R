# Model B: add_group_effect() -- nominal, post-hoc group effects.

make_wide <- function(n = 400, p = 60, seed = 1) {
  set.seed(seed)
  as.data.frame(matrix(rnorm(n * p), n, p,
                       dimnames = list(NULL, paste0("f", seq_len(p)))))
}

# ---- generative (spike-slab) mode ----------------------------------------

test_that("generative: adds factor with requested levels and proportions", {
  df  <- make_wide()
  out <- add_group_effect(df, name = "tissue",
                          levels = c("liver", "muscle", "fat"),
                          probs  = c(0.5, 0.3, 0.2),
                          effect_prop = 0.2, seed = 42)
  expect_true("tissue" %in% names(out))
  expect_s3_class(out$tissue, "factor")
  expect_equal(levels(out$tissue), c("liver", "muscle", "fat"))
  props <- as.numeric(prop.table(table(out$tissue)))
  expect_lt(max(abs(props - c(0.5, 0.3, 0.2))), 0.05)
})

test_that("generative: each level hits ~effect_prop of features, in DIFFERENT sets (nominal)", {
  df  <- make_wide(p = 200)
  out <- add_group_effect(df, name = "g", levels = c("a", "b", "c"),
                          effect_prop = 0.2, effect_sd = 0.5, seed = 7)
  sm  <- attr(out, "group_effect")$shift_matrix
  fracs <- rowMeans(sm != 0)
  expect_true(all(abs(fracs - 0.2) < 0.1))
  # nominal: the hit sets differ between levels (not the same columns)
  hit_a <- which(sm["a", ] != 0)
  hit_b <- which(sm["b", ] != 0)
  expect_false(identical(hit_a, hit_b))
})

test_that("generative: reference_level is left unshifted", {
  df  <- make_wide(p = 120)
  out <- add_group_effect(df, name = "g", levels = c("ctrl", "x", "y"),
                          effect_prop = 0.3, reference_level = "ctrl", seed = 3)
  sm  <- attr(out, "group_effect")$shift_matrix
  expect_true(all(sm["ctrl", ] == 0))
  expect_true(any(sm["x", ] != 0))
})

test_that("generative: integer target stays integer", {
  set.seed(9)
  df <- data.frame(count = as.integer(rpois(400, 20)),
                   z     = rnorm(400))
  out <- add_group_effect(df, name = "grp", levels = c("a", "b"),
                          on = "count", effect_prop = 1, effect_sd = 0.8, seed = 1)
  expect_true(is.integer(out$count))
})

# ---- explicit mode --------------------------------------------------------

test_that("explicit: per-level group means shift by the designed SD amount", {
  set.seed(2)
  df <- data.frame(x = rnorm(2000, 10, 4), y = rnorm(2000, 0, 2))
  out <- add_group_effect(df, name = "soil",
                          levels = c("clay", "sand", "silt"),
                          effects = list(clay = c(x = 1.0),    # +1 SD on x for clay
                                         silt = c(y = -0.8)),  # -0.8 SD on y for silt
                          seed = 5)
  sx <- sd(df$x); sy <- sd(df$y)
  # clay vs sand (reference) on x ~ +1 SD
  dx <- mean(out$x[out$soil == "clay"]) - mean(out$x[out$soil == "sand"])
  expect_lt(abs(dx / sx - 1.0), 0.15)
  # silt vs sand on y ~ -0.8 SD; x untouched for silt
  dy <- mean(out$y[out$soil == "silt"]) - mean(out$y[out$soil == "sand"])
  expect_lt(abs(dy / sy - (-0.8)), 0.15)
  # sand (reference) and omitted cells are zero in the matrix
  sm <- attr(out, "group_effect")$shift_matrix
  expect_equal(sm["sand", "x"], 0)
  expect_equal(sm["silt", "x"], 0)
})

test_that("explicit: bad level or feature names error", {
  df <- data.frame(x = rnorm(50), g = letters[1:50])
  expect_error(
    add_group_effect(df, name = "s", levels = c("a", "b"),
                     effects = list(zzz = c(x = 1))),
    "not in `levels`"
  )
  expect_error(
    add_group_effect(df, name = "s", levels = c("a", "b"),
                     effects = list(a = c(nope = 1))),
    "non-numeric or absent"
  )
})

# ---- shared behavior ------------------------------------------------------

test_that("works on a synthetica_sim and records $group_effect", {
  input <- data.frame(x = rnorm(300), y = rnorm(300), w = rnorm(300))
  sim   <- simulate_dataset(input, n = 300, seed = 42, verbose = FALSE)
  out   <- add_group_effect(sim, name = "soil", levels = c("clay", "sand"),
                            effect_prop = 0.5, seed = 11)
  expect_s3_class(out, "synthetica_sim")
  expect_true("soil" %in% names(out$data))
  expect_false(is.null(out$group_effect))
  expect_equal(out$group_effect$mode, "spike_slab")
  expect_equal(out$group_effect$levels, c("clay", "sand"))
})

test_that("group_col collision errors; non-target columns pass through", {
  df <- make_wide(p = 10)
  df$grp <- "x"
  expect_error(add_group_effect(df, name = "grp", levels = c("a", "b")),
               "collides")
  # character column passes through untouched
  df2 <- make_wide(p = 5); df2$lab <- sample(letters[1:3], nrow(df2), TRUE)
  out <- add_group_effect(df2, name = "g", levels = c("a", "b"),
                          effect_prop = 0.4, seed = 1)
  expect_equal(out$lab, df2$lab)
})

test_that("reproducible under seed", {
  df <- make_wide(p = 30)
  a  <- add_group_effect(df, name = "g", levels = c("a", "b", "c"),
                         effect_prop = 0.3, seed = 99)
  b  <- add_group_effect(df, name = "g", levels = c("a", "b", "c"),
                         effect_prop = 0.3, seed = 99)
  expect_equal(a, b)
})
