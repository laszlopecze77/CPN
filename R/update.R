#' Update a CPN Model with a New Formula or Data
#'
#' Updates a fitted Compound Poisson-Normal (CPN) model by changing the model
#' formula or data.
#'
#' @param object An object of class \code{"cpn"}, typically produced by a CPN
#' model fitting function.
#' @param new_formula A formula object specifying the new model formula.
#' @param data A data frame containing the new data. If \code{NULL},
#' the original data from the object is used.
#' @param ... Additional arguments passed to the underlying model fitting
#' function.
#'
#' @return A new \code{"cpn"} object fitted with the updated formula
#' and/or data.
#'
#' @seealso \code{\link{cpn}}
#'
#' @export
update.cpn <- function(object, new_formula, data = NULL, ...) {
  new_formula <- stats::update.formula(object$formula, new_formula)
  new_data <- if (!is.null(data)) data else object$data
  cpn(new_formula, data = new_data, ...)
}
