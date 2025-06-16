#' Bayesian Information Criterion for Compound Poisson-Normal Models
#'
#' Computes the Bayesian Information Criterion (BIC) for a fitted Compound
#' Poisson-Normal (CPN) regression model.
#'
#' @param object An object of class \code{"cpn"}, typically the result of a
#'   call to \code{\link{cpn}}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The BIC is computed as:
#' \deqn{-2 \cdot \log L + k \cdot \log(n)}
#' where \eqn{L} is the likelihood of the fitted model, \eqn{k} is the number
#' of estimated parameters (including regression coefficients, \eqn{\mu}, and
#' \eqn{\sigma}), and \eqn{n} is the number of observations.
#'
#' @return A numeric value representing the BIC of the fitted model.
#'
#' @seealso \code{\link{cpn}}, \code{\link{AIC}}, \code{\link[stats]{BIC}}
#'
#' @export
BIC.cpn <- function(object, ...) {            # nolint
  loglik <- -object$neg_log_likelihood
  n <- nrow(object$model)
  k <- length(object$se) # includes beta, mu, sigma
  bic <- -2 * loglik + k * log(n)
  bic
}
