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

library(testthat)

simulate_cpn_data <- function(n = 100, beta = c(0.5, -0.3, 0.7), mu = 1, sigma = 2) {
  set.seed(123)  # For reproducibility

  # Simulate predictors
  x1 <- factor(sample(c("A", "B"), size = n, replace = TRUE))  # Categorical
  x2 <- rnorm(n)  # Continuous

  # Create model matrix (includes intercept and dummy variable for x1)
  X <- model.matrix(~ x1 + x2)  # Will generate intercept, x1B, x2

  # Compute linear predictor and Poisson rates
  eta <- X %*% beta
  lambda <- exp(eta)

  # Simulate response variable
  y <- numeric(n)
  for (i in 1:n) {
    k <- rpois(1, lambda[i])
    if (k > 0) {
      y[i] <- sum(rnorm(k, mean = mu, sd = sigma))
    } else {
      y[i] <- 0
    }
  }

  # Return as a data.frame
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("cpn estimates are close to true parameters", {
  # Simulate data
  set.seed(123)
  true_beta <- c(0.5, -0.3, 0.7)
  true_mu <- 1
  true_sigma <- 2
  n <- 1000
  dat <- simulate_cpn_data(n = n, beta = true_beta, mu = true_mu, sigma = true_sigma)

  # Fit model
  fit <- cpn(y ~ x1 + x2, data = dat, k_max = 15)

  # Extract estimates
  beta_hat <- coef(fit, full=F)
  mu_hat <- fit$mu
  sigma_hat <- fit$sigma

  # Check beta names match expected
  expect_named(beta_hat, c("(Intercept)", "x1B", "x2"))

  # Compare estimated vs true
  expect_equal(as.numeric(beta_hat["(Intercept)"]), true_beta[1], tolerance = 0.1)
  expect_equal(as.numeric(beta_hat["x1B"]),        true_beta[2], tolerance = 0.1)
  expect_equal(as.numeric(beta_hat["x2"]),        true_beta[3], tolerance = 0.1)
  expect_equal(as.numeric(mu_hat),                true_mu,      tolerance = 0.1)
  expect_equal(as.numeric(sigma_hat),             true_sigma,   tolerance = 0.1)
})
