#' Predict Method for CPN Model Fits
#'
#' Generate predictions from a fitted Compound Poisson-Normal (CPN) regression model.
#'
#' @param object An object of class `"cpn"`, typically the result of a call to [cpn()].
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param type Type of prediction:
#'   \describe{
#'     \item{"link"}{Returns the linear predictor \eqn{\eta = X \beta}.}
#'     \item{"response"}{Returns the predicted mean response \eqn{E[Y] = \lambda \cdot \mu}, where \eqn{\lambda = \exp(\eta)} and \eqn{\mu} is the mean of the Normal component.}
#'   }
#' @param ... Currently ignored.
#'
#' @return A numeric vector of predicted values, either on the link scale (`type = "link"`) or the response scale (`type = "response"`).
#'
#' @seealso [cpn()] for model fitting.
#'
#' @examples
#' \dontrun{
#' fit <- cpn(y ~ x1 + x2, data = dat)
#' predict(fit, newdata = dat, type = "response")
#' predict(fit, newdata = dat, type = "link")
#' }
#'
#' @method predict cpn
#' @export

predict.cpn <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)

  # Remove response from terms to avoid 'y not found' error
  terms_noy <- stats::delete.response(object$terms)

  if (is.null(newdata)) {
    X <- stats::model.matrix(terms_noy, data = object$model)
  } else {
    mf <- stats::model.frame(terms_noy, data = newdata)
    X <- stats::model.matrix(terms_noy, data = mf)
  }

  beta_hat <- object$coefficients
  eta <- as.vector(X %*% beta_hat)

  if (type == "link") {
    return(eta)
  } else if (type == "response") {
    lambda_hat <- exp(eta)
    mu_hat <- object$mu
    return(lambda_hat * mu_hat)
  }
}


