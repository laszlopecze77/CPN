
test_that("vcov.cpn returns valid variance-covariance matrix", {

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
