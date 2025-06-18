#' Negative Log-Likelihood for Compound Poisson-Normal Regression
#'
#' Computes the negative log-likelihood for a regression model where
#' the response is assumed to follow a Compound Poisson-Normal distribution.
#'
#' @param beta_mu_sigma Numeric vector. Contains the regression coefficients
#'   for the log-link Poisson mean (`lambda`), followed by the mean (`mu`)
#'   and standard deviation (`sigma`) of the Normal components.
#' @param X Numeric matrix. The design matrix for the regression model
#'   (rows = observations, columns = covariates).
#' @param y Numeric vector. The observed responses.
#' @param k_max Integer. Maximum number of Poisson events to consider for
#'   truncation in likelihood approximation (default is 10).
#'
#' @return A numeric value giving the negative log-likelihood. Returns
#'   `Inf` if invalid parameters are provided or numerical issues occur.
#'
#' @details The function computes the likelihood for each observation as a
#'   weighted sum over a finite number of Poisson-Normal mixtures,
#'   truncating at `k_max`. If any resulting likelihood is zero, the function
#'   returns `Inf` to penalize the parameter set.
#'
#' @examples
#' X <- matrix(c(1, 2, 3, 4), ncol = 2)
#' y <- c(0.5, 1.5)
#' beta_mu_sigma <- c(0.1, 0.2, 1, 0.5)
#' cpn_neg_log_likelihood(beta_mu_sigma, X, y, k_max = 10)
#'
#' @export
cpn_neg_log_likelihood <- function(beta_mu_sigma,
                                   X, # nolint
                                   y,
                                   k_max = 10) {
  stopifnot(
    !is.na(beta_mu_sigma),
    is.numeric(beta_mu_sigma),
    is.numeric(X),
    is.numeric(y),
    is.numeric(k_max)
  )

  # Extract parameters
  p <- ncol(X)
  beta <- beta_mu_sigma[1:p]        # Regression coefficients for lambda
  mu <- beta_mu_sigma[p + 1]        # Mean of normal components
  sigma <- beta_mu_sigma[p + 2]     # Std dev of normal components

  if (sigma <= 0) return(Inf)       # Penalize invalid sigma

  # Linear predictor for lambda (log-link)
  eta <- X %*% beta
  lambda_vec <- as.vector(exp(eta))  # Observation-specific lambda

  # Compute likelihood for each observation
  compute_likelihood <- function(x_i, lambda_i) {
    if (x_i == 0) {
      stats::dpois(0, lambda_i)
    } else {
      k_vals <- seq_len(k_max)
      poisson_probs <- stats::dpois(k_vals, lambda_i)
      normal_probs <- stats::dnorm(
        x_i,
        mean = k_vals * mu,
        sd = sqrt(k_vals * sigma^2)
      )
      sum(poisson_probs * normal_probs)
    }
  }

  likelihoods <- mapply(compute_likelihood, x_i = y, lambda_i = lambda_vec)

  # If any likelihood is 0 (log(0) = -Inf), return Inf for the negative
  # log-likelihood
  if (any(likelihoods <= 0)) return(Inf)

  -sum(log(likelihoods))
}
