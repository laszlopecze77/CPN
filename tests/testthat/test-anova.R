library(testthat)

set.seed(123)
simulate_cpn_data <- function(n = 100,
                              beta = c(0.5, -0.3, 0.7),
                              mu = 1,
                              sigma = 2) {
  set.seed(123)  # For reproducibility

  # Simulate predictors
  x1 <- factor(sample(c("A", "B"), size = n, replace = TRUE))  # Categorical
  x2 <- rnorm(n)  # Continuous

  # Create model matrix (includes intercept and dummy variable for x1)
  X <- model.matrix(~ x1 + x2)  # Will generate intercept, x1B, x2  # nolint

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


test_that("anova.cpn works with a single model (Type I)", {

  fit <- cpn(y ~ x1 + x2, data = test_data)
  aov_tbl <- anova(fit)

  expect_s3_class(aov_tbl, "anova.cpn")
  expect_s3_class(aov_tbl, "data.frame")
  expect_true("Term" %in% names(aov_tbl))
  expect_true(any(!is.na(aov_tbl$`Pr(>Chi)`[-1])))  # should compute p-values
})

test_that("anova.cpn works with multiple nested models", {
  fit1 <- cpn(y ~ 1, data = test_data)
  fit2 <- cpn(y ~ x1, data = test_data)
  fit3 <- cpn(y ~ x1 + x2, data = test_data)


  comp_tbl <- anova(fit1, fit2, fit3)

  expect_s3_class(comp_tbl, "anova.cpn")
  expect_equal(nrow(comp_tbl), 3)
  expect_true("Model" %in% names(comp_tbl))
  expect_true("Pr(>Chi)" %in% names(comp_tbl))
  expect_true(all(is.na(comp_tbl$`Pr(>Chi)`[1])))  # First row should be NA
})

test_that("anova.cpn returns NULL for intercept-only model", {
  fit <- cpn(y ~ 1, data = test_data)

  expect_warning(out <- anova(fit), "only an intercept")
  expect_null(out)
})


test_that(
  "anova.cpn produces a valid sequential deviance table for a single model",
  {
    # Fit a full model with two predictors
    full_model <- cpn(y ~ x1 + x2, data = test_data)

    # Run ANOVA sequentially
    aov_seq <- anova(full_model)

    expect_s3_class(aov_seq, "anova.cpn")

    # There should be 1 row for residuals + 1 per term
    expect_equal(nrow(aov_seq), 3)
    expect_equal(aov_seq$Term[1], "Residuals")
    expect_equal(sort(aov_seq$Term[-1]), sort(c("x1", "x2")))

    # Residual deviance should decrease
    expect_true(
      all(diff(aov_seq$`Resid. Dev`[-1]) < 0)
    )

    # Residual degrees of freedom should decrease
    expect_true(
      all(diff(aov_seq$`Resid. Df`[-1]) < 0)
    )

    # Deviance differences should be non-negative
    expect_true(
      all(aov_seq$Deviance[-1] >= 0)
    )

    # Check that p-values are within [0, 1] or NA
    expect_true(
      all(
        is.na(aov_seq$`Pr(>Chi)`[-1]) |
          (aov_seq$`Pr(>Chi)`[-1] >= 0 & aov_seq$`Pr(>Chi)`[-1] <= 1)
      )
    )

    # First row (Residuals) should have NA for test columns
    expect_true(
      all(
        is.na(aov_seq[1, c("Df", "Deviance", "Pr(>Chi)")])
      )
    )
  }
)
