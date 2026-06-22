test_that("add_batch_effect adds a factor column with correct levels and balanced sizes", {
  set.seed(1)
  df  <- data.frame(x = rnorm(200), y = rnorm(200))
  out <- add_batch_effect(df, n_batches = 4, seed = 42)

  expect_true("batch" %in% names(out))
  expect_s3_class(out$batch, "factor")
  expect_equal(levels(out$batch), paste0("batch", 1:4))
  expect_equal(sum(table(out$batch)), 200L)
  # equally-sized to within 1
  expect_true(all(abs(table(out$batch) - 50L) <= 1L))
})

test_that("add_batch_effect: NULL/zero shifts leave numeric values unchanged", {
  set.seed(1)
  df  <- data.frame(x = rnorm(200), y = rnorm(200))
  out <- add_batch_effect(df, n_batches = 4, shifts = c(0, 0, 0, 0), seed = 42)
  expect_equal(out$x, df$x)
  expect_equal(out$y, df$y)
})

test_that("add_batch_effect: a large shift produces visible group separation", {
  set.seed(1)
  df  <- data.frame(x = rnorm(400), y = rnorm(400), z = rnorm(400))
  out <- add_batch_effect(df, n_batches = 4,
                          shifts = c(0, 0, 0, 2),
                          seed   = 42)
  # Batch 4 received a +2 SD shift; group means should differ by ~2 SD
  m1 <- mean(out$x[out$batch == "batch1"])
  m4 <- mean(out$x[out$batch == "batch4"])
  expect_gt(m4 - m1, 1.5)
})

test_that("add_batch_effect: works on a synthetica_sim object", {
  input <- data.frame(x = rnorm(200), y = rnorm(200), sex = sample(0:1, 200, TRUE))
  sim   <- simulate_dataset(input, n = 200, seed = 42, verbose = FALSE)
  out   <- add_batch_effect(sim, n_batches = 4, shifts = c(0, 0, 0, 1), seed = 7)

  expect_s3_class(out, "synthetica_sim")
  expect_true("batch" %in% names(out$data))
  expect_false(is.null(out$batch_effect))
  expect_equal(out$batch_effect$shifts, c(0, 0, 0, 1))
  expect_equal(out$batch_effect$n_batches, 4L)
})

test_that("add_batch_effect: batch_col name collision errors", {
  df <- data.frame(x = rnorm(50), batch = seq_len(50))
  expect_error(add_batch_effect(df, n_batches = 4), "collides")
})

test_that("add_batch_effect: non-numeric columns are passed through", {
  set.seed(1)
  df <- data.frame(
    x = rnorm(200),
    g = factor(sample(letters[1:3], 200, replace = TRUE))
  )
  out <- add_batch_effect(df, n_batches = 4, shifts = c(0, 0, 0, 1), seed = 42)
  expect_equal(out$g, df$g)
  expect_true("batch" %in% names(out))
})

test_that("add_batch_effect: custom level_names propagate", {
  set.seed(1)
  df  <- data.frame(x = rnorm(100))
  out <- add_batch_effect(df, n_batches = 3,
                          level_names = c("plate_A", "plate_B", "plate_C"),
                          seed = 1)
  expect_equal(levels(out$batch), c("plate_A", "plate_B", "plate_C"))
})

test_that("add_batch_effect: `on` accepts a character vector of feature names", {
  set.seed(1)
  df <- data.frame(x = rnorm(200), y = rnorm(200), z = rnorm(200))
  out <- add_batch_effect(df, n_batches = 4,
                          shifts = c(0, 0, 0, 2),
                          on     = c("x", "z"),
                          seed   = 42)
  # y should be unchanged; x and z should differ between batches
  expect_equal(out$y, df$y)
  expect_gt(mean(out$x[out$batch == "batch4"]) - mean(out$x[out$batch == "batch1"]), 1.5)
  expect_gt(mean(out$z[out$batch == "batch4"]) - mean(out$z[out$batch == "batch1"]), 1.5)
})

test_that("add_batch_effect: seed produces reproducible output", {
  df <- data.frame(x = rnorm(100), y = rnorm(100))
  a  <- add_batch_effect(df, n_batches = 4, shifts = c(0, 0, 0, 1), seed = 99)
  b  <- add_batch_effect(df, n_batches = 4, shifts = c(0, 0, 0, 1), seed = 99)
  expect_equal(a, b)
})

# ---- spike-and-slab generative model -------------------------------------

test_that("add_batch_effect: spike_slab perturbs roughly shift_prob of features", {
  set.seed(1)
  df  <- as.data.frame(matrix(rnorm(200 * 100), nrow = 200,
                              dimnames = list(NULL, paste0("f", 1:100))))
  out <- add_batch_effect(df, n_batches = 3, shifts = "spike_slab",
                          shift_prob = 0.2, shift_sd = 0.5, seed = 42)
  sm  <- attr(out, "batch_effect")$shift_matrix

  # reference batch (1) is fully unshifted
  expect_true(all(sm[1, ] == 0))

  # non-reference batches hit a fraction near shift_prob (absolute band)
  frac2 <- mean(sm[2, ] != 0)
  frac3 <- mean(sm[3, ] != 0)
  expect_lt(abs(frac2 - 0.2), 0.12)
  expect_lt(abs(frac3 - 0.2), 0.12)
})

test_that("add_batch_effect: reference_batch = NULL perturbs all batches", {
  set.seed(1)
  df  <- as.data.frame(matrix(rnorm(200 * 80), nrow = 200,
                              dimnames = list(NULL, paste0("f", 1:80))))
  out <- add_batch_effect(df, n_batches = 3, shifts = "spike_slab",
                          shift_prob = 0.3, reference_batch = NULL, seed = 7)
  sm  <- attr(out, "batch_effect")$shift_matrix
  # every batch should have at least one perturbed feature
  expect_true(all(rowSums(sm != 0) > 0))
})

test_that("add_batch_effect: shift_matrix is recorded with correct dimnames", {
  set.seed(1)
  df  <- data.frame(x = rnorm(200), y = rnorm(200), z = rnorm(200))
  out <- add_batch_effect(df, n_batches = 4, shifts = "spike_slab",
                          on = c("x", "z"), seed = 5)
  sm  <- attr(out, "batch_effect")$shift_matrix
  expect_equal(dim(sm), c(4L, 2L))
  expect_equal(rownames(sm), paste0("batch", 1:4))
  expect_equal(colnames(sm), c("x", "z"))
  # y was not a target, so it must be untouched
  expect_equal(out$y, df$y)
})

# ---- unequal batch sizes (proportions) -----------------------------------

test_that("add_batch_effect: proportions produce the requested batch sizes", {
  set.seed(1)
  df  <- data.frame(x = rnorm(1000))
  out <- add_batch_effect(df, n_batches = 4,
                          proportions = c(0.4, 0.3, 0.2, 0.1), seed = 42)
  tab <- as.integer(table(out$batch))
  expect_equal(sum(tab), 1000L)
  expect_lt(max(abs(tab - c(400, 300, 200, 100))), 2L)
})

test_that("add_batch_effect: bad proportions are rejected", {
  df <- data.frame(x = rnorm(100))
  expect_error(add_batch_effect(df, n_batches = 4, proportions = c(0.5, 0.5)),
               "length n_batches")
  expect_error(add_batch_effect(df, n_batches = 2, proportions = c(0.5, 0.6)),
               "sum to 1")
})

test_that("add_batch_effect: spike_slab is reproducible and stored on a sim", {
  input <- data.frame(x = rnorm(200), y = rnorm(200), w = rnorm(200))
  sim   <- simulate_dataset(input, n = 200, seed = 42, verbose = FALSE)
  a <- add_batch_effect(sim, n_batches = 3, shifts = "spike_slab", seed = 11)
  b <- add_batch_effect(sim, n_batches = 3, shifts = "spike_slab", seed = 11)
  expect_equal(a$data, b$data)
  expect_equal(a$batch_effect$model, "spike_slab")
  expect_equal(a$batch_effect$shift_prob, 0.2)
  expect_false(is.null(a$batch_effect$shift_matrix))
})
