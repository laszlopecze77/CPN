#' Extract coefficients from a Compound Poisson-Normal model
#'
#' Returns the estimated regression coefficients from a fitted Compound Poisson-Normal (CPN) model.
#' Optionally includes the estimated distribution parameters `mu` and `sigma`.
#'
#' @param object An object of class \code{cpn}.
#' @param full Logical. If \code{TRUE} (default), returns both the regression coefficients and the
#'   estimated distribution parameters \code{mu} and \code{sigma}. If \code{FALSE}, returns only the
#'   regression coefficients.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector of model coefficients. If \code{full = TRUE}, includes
#'   \code{mu} and \code{sigma}.
#' @export
#'
coef.cpn <- function(object, full = TRUE, ...) {
  if (full) {
    return(c(object$coefficients, object$mu, object$sigma))
  }
  return(object$coefficients)
}
