#' Variance-Covariance Matrix for a CPN Model
#'
#' Computes the variance-covariance matrix of parameter estimates from a fitted
#' Compound Poisson-Normal (CPN) regression model using the numerical Hessian
#' of the negative log-likelihood.
#'
#' @param object An object of class \code{"cpn"} returned
#' by the \code{\link{cpn}} function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A variance-covariance matrix of the model parameters. Rows and
#' columns are named
#' according to the model parameters (regression coefficients, \code{mu},
#' and \code{sigma}).
#' If the Hessian is singular, contains \code{NA} values, or is otherwise
#' invalid, a matrix
#' of \code{NA} values is returned with an appropriate warning.
#'
#' @details
#' The Hessian matrix of the negative log-likelihood is computed using
#' numerical finite differences
#' (via the \code{\link[numDeriv]{hessian}} function). The variance-covariance
#' matrix is then obtained
#' by inverting this Hessian. The result reflects local curvature and can be
#' used to compute standard
#' errors and confidence intervals for the parameters.
#'
#' @seealso \code{\link{cpn}}, \code{\link{summary.cpn}},
#'  \code{\link[numDeriv]{hessian}}
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(x = rnorm(100))
#' df$y <- sapply(exp(0.5 * df$x), function(lam) {
#'   k <- rpois(1, lam)
#'   if (k == 0) return(0)
#'   sum(rnorm(k, mean = 1, sd = 1))
#' })
#' fit <- cpn(y ~ x, data = df)
#' vcov(fit)
#'
#' @importFrom numDeriv hessian
#' @export
vcov.cpn <- function(object, ...) {
  # Full parameter vector (coefficients + mu + sigma)
  full_params <- c(object$coefficients, object$mu, object$sigma)

  # Compute Hessian
  H <- numDeriv::hessian(
    cpn_regression_neg_log_likelihood,
    x = full_params,
    X = stats::model.matrix(object$formula, object$data),
    y = stats::model.response(stats::model.frame(object$formula, object$data))
  )

  # Check Hessian validity
  if (any(is.na(H)) || det(H) == 0 || any(!is.finite(H))) {
    warning("Hessian is singular or contains non-finite values; vcov
            not available.")
    return(matrix(NA, nrow = length(full_params), ncol = length(full_params),
                  dimnames = list(names(full_params), names(full_params))))
  }

  # Return named variance-covariance matrix
  V <- solve(H)
  dimnames(V) <- list(names(full_params), names(full_params))
  return(V)
}
