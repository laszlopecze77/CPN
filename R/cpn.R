#' Compound Poisson-Normal Regression
#'
#' Fits a Compound Poisson-Normal (CPN) regression model to the response
#' variable using maximum likelihood estimation via the Nelder-Mead method
#' \insertCite{Nelder1965}{CPN}.
#'
#' @param formula A formula specifying the model, e.g., `y ~ x1 + x2`.
#' @param data An optional data frame, list, or environment containing the
#'   variables in the model. If not found in `data`, the variables are
#'   taken from the environment from which `cpn` is called.
#' @param mu_init Optional initial value for the Normal mean parameter.
#'   If `NULL`, it is set to the sample mean of `y`.
#' @param sigma_init Optional initial value for the Normal standard deviation.
#'   If `NULL`, it is set to the sample standard deviation of `y`.
#' @param k_max Upper limit of summation used in approximating the Poisson
#'   convolution. Should be between 10 and 100. Default is 10.
#'
#' @details
#' The function fits a regression model where the response is assumed to follow
#' a Compound Poisson-Normal distribution. The model estimates regression
#' coefficients for the Poisson intensity, and mean and standard deviation of
#' the Normal component.
#'
#' The optimization is performed via maximum likelihood using the Nelder-Mead
#' method. Standard errors are computed via the observed information matrix
#' using numerical Hessian. If the Hessian is singular or not positive
#' definite, a warning is issued and standard errors are returned as `NA`.
#'
#' @return An object of class `"cpn"` containing the following components:
#' \item{coefficients}{Estimated regression coefficients.}
#' \item{mu}{Estimated mean of the Normal component.}
#' \item{sigma}{Estimated standard deviation of the Normal component.}
#' \item{se}{Standard errors of the estimated parameters.}
#' \item{model}{The model frame used.}
#' \item{terms}{The terms object from the model.}
#' \item{formula}{The model formula.}
#' \item{data}{The original data passed in.}
#' \item{deviance_residuals}{Vector of deviance residuals.}
#' \item{fitted_values}{Fitted values (Poisson mean times Normal mean).}
#' \item{neg_log_likelihood}{Negative log-likelihood at the optimum.}
#' \item{null_deviance}{Null model deviance.}
#' \item{residual_deviance}{Residual deviance of the fitted model.}
#' \item{df_null}{Degrees of freedom for the null model.}
#' \item{df_residual}{Degrees of freedom for the fitted model.}
#' \item{k_max}{Value used to truncate the Poisson convolution sum.}
#' \item{call}{The matched function call.}
#'
#' @references
#' \insertRef{Nelder1965}{CPN}
#'
#' @examples
#' set.seed(123)
#' data <- simulate_cpn_data()
#'
#' # Sequential analysis of deviance
#' fit <- cpn(y ~ x1 + x2, data = data)
#' summary(fit)
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats as.formula model.frame model.matrix model.response
#' @importFrom stats optim setNames sd dnorm dpois terms
#' @export
#'
cpn <- function(formula,
                data = NULL,
                mu_init = NULL,
                sigma_init = NULL,
                k_max = 10) {

  null_or_numeric_1l <- function(x) {
    is.null(x) | (is.numeric(x) & length(x) == 1L)
  }
  stopifnot(
    inherits(formula, "formula"),
    null_or_numeric_1l(mu_init),
    null_or_numeric_1l(sigma_init),
    is.null(sigma_init) | sigma_init >= 0,
    null_or_numeric_1l(k_max)
  )
  if (k_max < 10 || k_max > 100) stop("k_max should be between 10 and 100")

  # Build model frame and design matrix
  # (allow formula to work with or without data argument)
  formula <- stats::as.formula(formula, env = parent.frame())
  mf <- stats::model.frame(formula = formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(attr(mf, "terms"), data = mf) # nolint

  # Initial values
  if (is.null(mu_init)) mu_init <- mean(y)
  if (is.null(sigma_init)) sigma_init <- stats::sd(y)
  init_beta <- rep(0, ncol(X))
  init_vals <- c(init_beta, mu_init, sigma_init)

  # Optimize negative log-likelihood
  fit <- stats::optim(
    par = init_vals,
    fn = cpn_neg_log_likelihood,
    X = X,
    y = y,
    k_max = k_max,
    method = "Nelder-Mead",
    control = list(maxit = 1000)
  )

  beta_hat <- fit$par
  loglik <- -fit$value

  # Compute Hessian and standard errors

  H <- tryCatch({ # nolint
    numDeriv::hessian(
      func = function(par) {
        cpn_neg_log_likelihood(
          par, X, y, k_max = k_max
        )
      },
      x = beta_hat
    )
  }, error = function(e) {
    warning("Failed to compute Hessian: ", e$message)
    matrix(NA, length(beta_hat), length(beta_hat))
  })

  if (any(is.na(H)) || det(H) == 0 || any(!is.finite(H))) {
    warning("Hessian is singular or contains non-finite values; SEs are not
            available.")
    se_hat <- rep(NA, length(beta_hat))
  } else {
    eigs <- eigen(H, symmetric = TRUE)$values
    if (any(eigs <= 0)) {
      warning("Hessian is not positive definite; SEs may be unreliable.")
      se_hat <- sqrt(diag(solve(H)))
    } else {
      se_hat <- sqrt(diag(solve(H)))
    }
  }

  param_names <- c(colnames(X), "mu", "sigma")
  names(beta_hat) <- param_names
  names(se_hat) <- param_names

  # Compute fitted values (linear predictor)
  eta <- as.vector(X %*% beta_hat[seq_len(ncol(X))])
  lambda_hat <- exp(eta) # Poisson mean
  mu_hat <- beta_hat["mu"]
  sigma_hat <- beta_hat["sigma"]

  # Compute individual log-likelihood contributions
  loglik_obs <- function(y_i,
                         lambda_i,
                         mu_hat,
                         sigma_hat,
                         k_max = k_max) {
    if (y_i == 0) {
      log(dpois(0, lambda_i))
    } else {
      k_vals <- 1:k_max
      pois_weights <- dpois(k_vals, lambda_i)
      norm_vals <- dnorm(
        y_i,
        mean = k_vals * mu_hat,
        sd = sqrt(k_vals) * sigma_hat
      )
      log(sum(pois_weights * norm_vals))
    }
  }


  dev_res <- numeric(length(y))
  for (i in seq_along(y)) {
    ll_hat <- loglik_obs(y[i], lambda_hat[i], mu_hat, sigma_hat, k_max)
    ll_sat <- 0
    diff <- y[i] - lambda_hat[i] * mu_hat
    dev_res[i] <- sign(diff) * sqrt(2 * (ll_sat - ll_hat))
  }

  # Null model (intercept only)
  X_null <- matrix(1, nrow = nrow(X), ncol = 1) # nolint
  colnames(X_null) <- "(Intercept)" # nolint
  null_fit <- stats::optim(
    par = c(0, mu_init, sigma_init),
    fn = cpn_neg_log_likelihood,
    X = X_null,
    y = y,
    method = "Nelder-Mead",
    k_max = k_max,
    control = list(maxit = 1000)
  )
  null_loglik <- -null_fit$value

  null_deviance <- -2 * null_loglik
  residual_deviance <- -2 * loglik
  df_null <- length(y) - 1
  df_residual <- length(y) - length(beta_hat)

  model_frame <- stats::model.frame(formula, data)
  terms_obj <- stats::terms(model_frame)

  structure(
    list(
      coefficients = stats::setNames(beta_hat[seq_len(ncol(X))], colnames(X)),
      mu = beta_hat["mu"],
      sigma = beta_hat["sigma"],
      se = se_hat,
      model = model_frame,
      terms = terms_obj,
      formula = formula,
      data = data,
      deviance_residuals = dev_res,
      fitted_values = lambda_hat * mu_hat,
      neg_log_likelihood = fit$value,
      null_deviance = null_deviance,
      residual_deviance = residual_deviance,
      df_null = df_null,
      df_residual = df_residual,
      k_max = k_max,
      call = match.call()
    ),
    class = "cpn"
  )
}
