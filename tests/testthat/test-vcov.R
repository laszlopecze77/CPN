
test_that("vcov.cpn returns valid variance-covariance matrix", {

  # Simulate test data
  simulate_cpn_data <- function(n = 100,
                                beta = c(0.5, -0.3, 0.7),
                                mu = 1,
                                sigma = 2) {
    # Simulate predictors
    x1 <- factor(sample(c("A", "B"), size = n, replace = TRUE))  # Categorical
    x2 <- rnorm(n)  # Continuous

    # Create model matrix (includes intercept and dummy variable for x1)
    X <- model.matrix(~ x1 + x2)  # Will generate intercept, x1B, x2 nolint

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

  set.seed(123)
  test_data <- simulate_cpn_data()

  # Fit the model
  fit <- cpn(y ~ x1 + x2, data = test_data)

  # Extract variance-covariance matrix
  V <- vcov(fit) # nolint

  # Check it's a matrix
  expect_true(is.matrix(V))

  # Check dimensions: should be (number of parameters) x (number of parameters)
  p <- length(fit$se)
  expect_equal(dim(V), c(p, p))

  # Check row and column names
  expect_equal(rownames(V), names(fit$se))
  expect_equal(colnames(V), names(fit$se))

  # Check that diagonal entries (variances) are non-negative
  expect_true(all(diag(V) >= 0, na.rm = TRUE))

  # Check that it's symmetric
  expect_equal(V, t(V))
})
