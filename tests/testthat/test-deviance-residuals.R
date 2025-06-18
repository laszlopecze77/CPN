# Assumes cpn() is already defined in your environment

set.seed(123)
n <- 100
x <- rnorm(n)
lambda <- exp(0.1 + 0.4 * x)
mu_true <- 1
sigma_true <- 2


k <- rpois(n, lambda)
y <- numeric(n)



for (i in 1:n) {
  if (k[i] == 0) {
    y[i] <- 0
  } else {
    y[i] <- sum(rnorm(k[i], mu_true, sigma_true))
  }
}

test_data <- data.frame(y = y, x = x)

# Fit the model
fit <- cpn(y ~ x, data = test_data)

# === TESTS ===

# 1. Residual deviance equals sum of squared deviance residuals
test_that("Residual deviance matches sum of squared deviance residuals", {
  expect_equal(sum(fit$deviance_residuals^2),
               fit$residual_deviance,
               tolerance = 1e-6)
})


# 2. Signs of residuals are consistent with observed - fitted
test_that(
  "Signs of residuals are consistent with observed - fitted", {
    res_signs <- sign(y - fit$fitted_values)
    dev_signs <- sign(fit$deviance_residuals)
    agreement <- mean(res_signs == dev_signs)
    expect_true(agreement > 0.98)
  }
)
