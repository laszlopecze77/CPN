#' Extract Log-Likelihood from a CPN Model
#'
#' Returns the log-likelihood of a fitted Compound Poisson-Normal (CPN) model.
#'
#' @param object An object of class \code{"cpn"}, typically resulting from a CPN model fitting function.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"logLik"} representing the log-likelihood.
#'   The \code{df} attribute contains the number of estimated parameters.
#'
#' @seealso \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{logLik}}
#'
#' @export
logLik.cpn <- function(object, ...) {
  val <- -object$neg_log_likelihood
  attr(val, "df") <- length(object$se)
  class(val) <- "logLik"
  return(val)
}
