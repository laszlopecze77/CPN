#' Plot Diagnostics for a CPN Model
#'
#' Generates diagnostic plots for a fitted Compound Poisson-Normal (CPN) model
#' object.
#' Options include a residuals vs fitted values plot and a Q-Q plot of the
#' deviance residuals.
#'
#' @param x An object of class \code{"cpn"}, typically resulting from a CPN
#' model fitting function.
#' @param which A character string specifying the type of plot to produce.
#'   Options are \code{"residuals"} (default) for a residuals vs fitted values
#'   plot,
#'   and \code{"qq"} for a Q-Q plot of deviance residuals.
#' @param ... Additional graphical parameters passed to the underlying plotting
#' functions.
#'
#' @details
#' The residuals vs fitted plot helps assess non-linearity, unequal error
#' variances, and outliers.
#' The Q-Q plot checks for normality of deviance residuals.
#'
#' @return This function is called for its side effects and does not return a
#' value.
#'
#' @seealso \code{\link{residuals.cpn}}, \code{\link{fitted.cpn}}
#'
#' @export

plot.cpn <- function(x, which = c("residuals", "qq"), ...) {
  which <- match.arg(which)
  graphics::par(mfrow = c(1, 1))
  if (which == "residuals") {
    plot(x$fitted_values, x$deviance_residuals,
         xlab = "Fitted values",
         ylab = "Deviance residuals",
         main = "Residuals vs Fitted")
    graphics::abline(h = 0, col = "red", lty = 2)
  } else if (which == "qq") {
    stats::qqnorm(x$deviance_residuals)
    stats::qqline(x$deviance_residuals, col = "red")
  }
}
