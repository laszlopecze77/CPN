if (requireNamespace("lintr", quietly = TRUE)) {
  test_that("Package is lint-free", {
    lintr::expect_lint_free()
  })
}
