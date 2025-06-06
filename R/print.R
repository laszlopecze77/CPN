#' Print method for Compound Poisson-Normal model objects
#'
#' Displays a concise summary of a fitted Compound Poisson-Normal (CPN)
#' regression model, including
#' the original function call, estimated coefficients, distribution
#' parameters (`mu` and `sigma`),
#' residual deviance, degrees of freedom, and AIC.
#'
#' @param x An object of class \code{cpn}.
#' @param ... Additional arguments passed to other methods (currently ignored).
#'
#' @return The input object \code{x}, invisibly.
#' @importFrom stats AIC
#' @export
#'

print.cpn <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  print(round(x$coefficients, 4))

  cat(sprintf("\nmu: %.4f\n", x$mu))
  cat(sprintf("sigma: %.4f\n", x$sigma))

  cat(sprintf("\nResidual deviance: %.2f on %d degrees of freedom\n",
              x$residual_deviance, x$df_residual))
  cat(sprintf("AIC: %.2f\n", AIC(x)))

  invisible(x)
}
