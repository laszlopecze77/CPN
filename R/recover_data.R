
#' Recover Data for `cpn` Model Objects
#'
#' This method is used by the \pkg{emmeans} package to recover the data and
#' terms from a fitted Compound Poisson-Normal (CPN) model object. It is
#' typically used internally when computing estimated marginal means.
#'
#' @param object A fitted model object of class `cpn`.
#' @param ... Additional arguments (not used).
#'
#' @return An object containing the recovered data and terms, suitable for use
#'   by \pkg{emmeans}.
#'
#' @importFrom emmeans recover_data
#' @export
#' @method recover_data cpn
#'
recover_data.cpn <- function(object, ...) {
  emmeans::recover_data(
    object$call,
    trms = stats::delete.response(object$terms),
    na.action = stats::na.omit,
    data = object$data
  )
}
