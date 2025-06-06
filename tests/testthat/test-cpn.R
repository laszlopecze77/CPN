test_that("cpn checks mu_init", {
  set.seed(123)
  n <- 3
  x <- rnorm(n)
  lambda <- exp(0.5 + 0.2 * x)
  N <- rpois(n, lambda) # nolint
  y <- ifelse(N == 0, 0, rnorm(n, mean = N * 3, sd = sqrt(N) * 2))

  expect_error(cpn(y ~ x, mu_init = c(1, 2)))
  expect_error(cpn(y ~ x, mu_init = "a"))
  expect_no_error(cpn(y ~ x, mu_init = 1))
  expect_no_error(cpn(y ~ x, mu_init = NULL))
  expect_no_error(cpn(y ~ x))
})

test_that("cpn checks sigma_init", {
  set.seed(123)
  n <- 3
  x <- rnorm(n)
  lambda <- exp(0.5 + 0.2 * x)
  N <- rpois(n, lambda) # nolint
  y <- ifelse(N == 0, 0, rnorm(n, mean = N * 3, sd = sqrt(N) * 2))

  expect_error(cpn(y ~ x, sigma_init = c(1, 2)))
  expect_error(cpn(y ~ x, sigma_init = "a"))
  expect_error(cpn(y ~ x, sigma_init = -1))

  expect_no_error(cpn(y ~ x, sigma_init = 1))
  expect_no_error(cpn(y ~ x, sigma_init = NULL))
  expect_no_error(cpn(y ~ x))
})

test_that("cpn checks k_max", {
  set.seed(123)
  n <- 3
  x <- rnorm(n)
  lambda <- exp(0.5 + 0.2 * x)
  N <- rpois(n, lambda) # nolint
  y <- ifelse(N == 0, 0, rnorm(n, mean = N * 3, sd = sqrt(N) * 2))

  expect_error(cpn(y ~ x, k_max = c(1, 2)))
  expect_error(cpn(y ~ x, k_max = 1))
  expect_error(cpn(y ~ x, k_max = 1000))
  expect_error(cpn(y ~ x, k_max = "a"))

  expect_no_error(cpn(y ~ x, k_max = 10))
  expect_no_error(cpn(y ~ x, k_max = 100))
  expect_no_error(cpn(y ~ x))
})

test_that("cpn fetch formula terms in env when data is not provided", {
  set.seed(123)
  n <- 3
  x <- rnorm(n)
  lambda <- exp(0.5 + 0.2 * x)
  N <- rpois(n, lambda) # nolint
  y <- ifelse(N == 0, 0, rnorm(n, mean = N * 3, sd = sqrt(N) * 2))
  expect_no_error(cpn(y ~ x, mu_init = 1))
  expect_error(cpn(y ~ x2, mu_init = 1))
})
