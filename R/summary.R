#' Summarize a Compound Poisson-Normal (CPN) Model Fit
#'
#' Produces a summary for a fitted `cpn` model object, including parameter
#' estimates,
#' standard errors, z-values, p-values, deviance residual summaries, and model
#' diagnostics.
#'
#' @param object An object of class `"cpn"`, typically resulting from a call
#' to a CPN regression function.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `"summary.cpn"`, which is a list containing:
#' \describe{
#'   \item{call}{The matched call that generated the model.}
#'   \item{summary_table}{A data frame with parameter estimates, standard
#'   errors, z-values, and p-values.}
#'   \item{deviance_summary}{A five-number summary (Min, 1Q, Median, 3Q, Max)
#'   of deviance residuals.}
#'   \item{mu}{Estimated value of the Normal component mean.}
#'   \item{sigma}{Estimated value of the Normal component standard deviation.}
#'   \item{null_deviance}{Null deviance of the model.}
#'   \item{residual_deviance}{Residual deviance of the model.}
#'   \item{df_null}{Degrees of freedom for the null model.}
#'   \item{df_residual}{Degrees of freedom for the residual model.}
#'   \item{aic}{Akaike Information Criterion for model comparison.}
#' }
#'
#' @seealso [cpn()], [anova.cpn()]
#' @method summary cpn
#' @export
summary.cpn <- function(object, ...) {
  beta_hat <- object$coefficients
  se_hat <- object$se[seq_along(beta_hat)]

  z_vals <- beta_hat / se_hat
  p_vals <- 2 * stats::pnorm(-abs(z_vals))
  param_names <- names(beta_hat)

  summary_table <- data.frame(
    Estimate = beta_hat,
    Std.Error = se_hat,
    z.value = z_vals,
    Pr.z = p_vals,
    row.names = param_names
  )

  dev_res <- object$deviance_residuals
  dev_summary <- stats::quantile(dev_res, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(dev_summary) <- c("Min", "1Q", "Median", "3Q", "Max")

  out <- list(
    call = object$call,
    summary_table = summary_table,
    deviance_summary = dev_summary,
    mu = object$mu,
    sigma = object$sigma,
    null_deviance = object$null_deviance,
    residual_deviance = object$residual_deviance,
    df_null = object$df_null,
    df_residual = object$df_residual,
    aic = object$aic
  )
  class(out) <- "summary.cpn"
  out
}

#' Print Method for Summary of CPN Model
#'
#' Displays a formatted summary of a fitted `cpn` model object, including
#' parameter estimates,
#' deviance residuals, and model fit statistics.
#'
#' @param x An object of class `"summary.cpn"`, as returned by [summary.cpn()].
#' @param ... Additional arguments passed to `print`, ignored by default.
#'
#' @return Invisibly returns the `summary.cpn` object.
#' @method print summary.cpn
#' @export
print.summary.cpn <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nDeviance Residuals:\n")
  print(x$deviance_summary, digits = 4, quote = FALSE)

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$summary_table, P.values = TRUE, has.Pvalue = TRUE)

  cat(sprintf("\nEstimated mu parameter: %.4f\n", x$mu))
  cat(sprintf("Estimated sigma parameter: %.4f\n", x$sigma))

  cat(sprintf("\nNull deviance: %.2f on %d degrees of freedom\n",
              x$null_deviance, x$df_null))
  cat(sprintf("Residual deviance: %.2f on %d degrees of freedom\n",
              x$residual_deviance, x$df_residual))
  cat(sprintf("AIC: %.2f\n", x$aic))

  invisible(x)
}
