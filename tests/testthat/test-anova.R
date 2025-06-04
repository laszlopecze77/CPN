library(testthat)


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


test_that("anova.cpn works with a single model (Type I)", {

  fit <- cpn(y ~ x, data = test_data)
  aov_tbl <- anova(fit)

  expect_s3_class(aov_tbl, "anova.cpn")
  expect_s3_class(aov_tbl, "data.frame")
  expect_true("Term" %in% names(aov_tbl))
  expect_true(any(!is.na(aov_tbl$`Pr(>Chi)`[-1])))  # should compute p-values
})

test_that("anova.cpn works with multiple nested models", {
  fit1 <- cpn(y ~ 1, data = test_data)
  fit2 <- cpn(y ~ x, data = test_data)


  comp_tbl <- anova(fit1, fit2)

  expect_s3_class(comp_tbl, "anova.cpn")
  expect_equal(nrow(comp_tbl), 2)
  expect_true("Model" %in% names(comp_tbl))
  expect_true("Pr(>Chi)" %in% names(comp_tbl))
  expect_true(all(is.na(comp_tbl$`Pr(>Chi)`[1])))  # First row should be NA
})

test_that("anova.cpn returns NULL for intercept-only model", {
  fit <- cpn(y ~ 1, data = test_data)

  expect_warning(out <- anova(fit), "only an intercept")
  expect_null(out)
})

