#' Predict Method for CPN Model Objects
#'
#' Computes predictions from a fitted Compound Poisson-Normal (CPN) regression
#' model.
#' Supports predictions on the link or response scale, with optional confidence
#' intervals.
#'
#' @param object An object of class `cpn`, typically the result of a call to a
#' function fitting a Compound Poisson-Normal regression model.
#' @param newdata An optional data frame in which to look for variables
#' with which to predict. If omitted, the original model data is used.
#' @param type Type of prediction: \code{"link"} returns the linear
#' predictor \eqn{\eta = X\beta}; \code{"response"} returns the mean
#' response \eqn{E[Y] = \mu \cdot \exp(\eta)}.
#' @param interval Type of interval calculation. Either \code{"none"} (default)
#' or \code{"confidence"} for confidence intervals around the predicted values.
#' @param level Confidence level for the interval. Defaults to 0.95.
#' @param ... Further arguments passed to or from other methods (not currently
#' used).
#'
#' @return If \code{interval = "none"}, returns a numeric vector of predicted
#' values on the specified scale.
#' If \code{interval = "confidence"}, returns a data frame with columns:
#' \describe{
#'   \item{\code{fit}}{Predicted value}
#'   \item{\code{lwr}}{Lower bound of the confidence interval}
#'   \item{\code{upr}}{Upper bound of the confidence interval}
#' }
#'
#' @details For predictions on the response scale with confidence intervals,
#' the standard errors of both the linear predictor and the estimated \code{mu}
#' parameter are combined using the delta method.
#'
#' Factor levels in \code{newdata} are aligned to match those used in the
#' original model fit.
#'
#' @seealso \code{\link{cpn}}, \code{\link{vcov}}, \code{\link{model.matrix}},
#'  \code{\link{predict}}
#'
#' @examples
#' \dontrun{
#' fit <- cpn(y ~ x1 + x2, data = mydata)
#' predict(fit, newdata = mynewdata, type = "response", interval = "confidence")
#' }
#'
#' @export

predict.cpn <- function(object, newdata = NULL, type = c("response", "link"),
                        interval = c("none", "confidence"), level = 0.95, ...) {


  type <- match.arg(type)
  interval <- match.arg(interval)

  terms_noy <- stats::delete.response(object$terms)

  if (is.null(newdata)) {
    X <- stats::model.matrix(terms_noy, data = object$model,     # nolint
                             contrasts.arg = attr(object$model, "contrasts"))
  } else {
    for (v in names(object$model)) {
      if (is.factor(object$model[[v]]) && v %in% names(newdata)) {
        newdata[[v]] <- factor(newdata[[v]], levels = levels(object$model[[v]]))
      }
    }
    mf <- stats::model.frame(terms_noy, data = newdata,
                             xlev = stats::.getXlevels(object$terms,
                                                       object$model))
    X <- stats::model.matrix(terms_noy, data = mf,            # nolint
                             contrasts.arg = attr(object$model, "contrasts"))
  }

  beta_hat <- object$coefficients
  eta <- as.vector(X %*% beta_hat)
  mu_hat <- object$mu
  se_mu <- object$se["mu"]
  z <- stats::qnorm(1 - (1 - level) / 2)

  if (interval == "none") {
    if (type == "link") return(eta)
    if (type == "response") return(exp(eta) * mu_hat)
  }

  # Confidence intervals
  V_full <- stats::vcov(object)     # nolint
  beta_names <- names(object$coefficients)  # only beta, not mu/sigma
  V_beta <- V_full[beta_names, beta_names, drop = FALSE] # nolint
  se_eta <- sqrt(rowSums((X %*% V_beta) * X))

  if (type == "link") {
    lwr <- eta - z * se_eta
    upr <- eta + z * se_eta
    return(data.frame(fit = eta, lwr = lwr, upr = upr))
  }

  # For response: include SE(mu) in delta method
  exp_eta <- exp(eta)
  fit <- exp_eta * mu_hat
  se_fit <- sqrt(
    (mu_hat * exp_eta)^2 * se_eta^2 +
      (exp_eta)^2 * se_mu^2
  )

  lwr <- fit - z * se_fit
  upr <- fit + z * se_fit
  return(data.frame(fit = fit, lwr = lwr, upr = upr))
}
