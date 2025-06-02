#' Extract Residuals from a Compound Poisson-Normal (CPN) Model
#'
#' This method extracts residuals from a fitted model object of
#' class \code{"cpn"}.
#' Supports both raw and deviance residuals.
#'
#' @param object An object of class \code{"cpn"}, typically produced by a
#' CPN model fitting function.
#' @param type Character string indicating the type of residuals to return.
#'   Choices are \code{"deviance"} (default) or \code{"raw"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of residuals. Returns deviance
#' residuals if \code{type = "deviance"},
#'   or raw residuals (observed - fitted) if \code{type = "raw"}.
#'
#' @export
residuals.cpn <- function(object, type = c("deviance", "raw"), ...) {
  type <- match.arg(type)

  if (type == "deviance") {
    return(object$deviance_residuals)
  } else if (type == "raw") {
    response <- object$data[[as.character(object$formula[[2]])]]
    return(response - object$fitted_values)
  }
}
