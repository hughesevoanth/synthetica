# Small synthetic input shared across tests. Pure simulation -- no real data.
#
# Default fixture is numeric-only so phenotype-injection tests using `on = "all"` don't emit
# the "dropped categorical target" warning. Pass `with_categorical = TRUE` for
# tests that specifically need a factor column (e.g., level-collapse tests).
make_test_input <- function(n = 200, seed = 1, with_categorical = FALSE) {
  set.seed(seed)
  d <- data.frame(
    x   = rnorm(n, 50, 10),
    y   = rlnorm(n, log(20), 0.5),
    sex = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.5, 0.5)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(with_categorical)) {
    d$grp <- factor(sample(letters[1:3], n, replace = TRUE))
  }
  d
}

# Numeric-only fixture for clean row-hashing in privacy tests.
make_numeric_fixture <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(
    x   = rnorm(n, 50, 10),
    y   = rlnorm(n, log(20), 0.5),
    sex = sample(c(0L, 1L), n, replace = TRUE)
  )
}
