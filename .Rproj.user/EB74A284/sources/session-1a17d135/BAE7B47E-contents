#' Compute Akaike Information Criterion (AIC) for a CPN Model
#'
#' Calculates the AIC for a fitted Compound Poisson-Normal (CPN) model.
#'
#' @param object An object of class \code{"cpn"}, typically produced by a CPN model fitting function.
#' @param ... Additional arguments (currently unused).
#' @param k Numeric penalty per parameter; the default is \code{2}, corresponding to the traditional AIC.
#'
#' @return A numeric value representing the AIC of the model.
#'
#' @seealso \code{\link{logLik}}, \code{\link{BIC}}, \code{\link{AIC}}
#'
#' @export

AIC.cpn <- function(object, ..., k = 2) {
  k * length(object$se) + 2 * object$neg_log_likelihood
}
