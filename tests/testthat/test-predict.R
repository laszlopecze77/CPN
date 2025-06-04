library(testthat)

set.seed(123)
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


test_data <- simulate_cpn_data()


test_that("predict.cpn returns expected values without newdata", {

  fit <- cpn(y ~ x1 + x2, data = test_data)

  # Response, no interval
  pred_resp <- predict(fit, type = "response")
  expect_type(pred_resp, "double")
  expect_length(pred_resp, nrow(test_data))

  # Link, no interval
  pred_link <- predict(fit, type = "link")
  expect_type(pred_link, "double")
  expect_length(pred_link, nrow(test_data))

  # Response with confidence interval
  pred_resp_ci <- predict(fit, type = "response", interval = "confidence")
  expect_s3_class(pred_resp_ci, "data.frame")
  expect_equal(names(pred_resp_ci), c("fit", "lwr", "upr"))
  expect_equal(nrow(pred_resp_ci), nrow(test_data))
  expect_true(all(pred_resp_ci$lwr <= pred_resp_ci$fit))
  expect_true(all(pred_resp_ci$upr >= pred_resp_ci$fit))

  # Link with confidence interval
  pred_link_ci <- predict(fit, type = "link", interval = "confidence")
  expect_s3_class(pred_link_ci, "data.frame")
  expect_equal(names(pred_link_ci), c("fit", "lwr", "upr"))
  expect_equal(nrow(pred_link_ci), nrow(test_data))
})

test_that("predict.cpn works with newdata", {
  set.seed(123)
  test_data <- simulate_cpn_data()
  fit <- cpn(y ~ x1 + x2, data = test_data)

  newdata <- test_data[1:10, ]
  pred_new <- predict(fit, newdata = newdata, type = "response")
  expect_type(pred_new, "double")
  expect_length(pred_new, 10)

  pred_new_ci <- predict(fit, newdata = newdata,
                         type = "response", interval = "confidence", level = 0.99)
  expect_s3_class(pred_new_ci, "data.frame")
  expect_equal(nrow(pred_new_ci), 10)
  expect_true(all(pred_new_ci$lwr <= pred_new_ci$fit))
  expect_true(all(pred_new_ci$upr >= pred_new_ci$fit))
})

test_that("predict.cpn handles factor levels correctly in newdata", {
  set.seed(123)
  test_data <- simulate_cpn_data()
  fit <- cpn(y ~ x1 + x2, data = test_data)

  newdata <- test_data[1:5, ]
  newdata$x1 <- factor(newdata$x1, levels = rev(levels(test_data$x1)))  # reversed order
  expect_silent(pred <- predict(fit, newdata = newdata))
  expect_length(pred, 5)
})
