#' Basis for Estimated Marginal Means for Compound Poisson-Normal Models
#'
#' Provides the required components for computing estimated marginal means
#' (EMMs) from a fitted Compound Poisson-Normal (CPN) model. This method is
#' used by the \pkg{emmeans} package to interface with models of class `"cpn"`.
#'
#' @param object A fitted model object of class `"cpn"` returned by a call to
#'   e.g., `cpn5()`.
#' @param trms Terms object extracted from the model formula, typically
#'   provided by \pkg{emmeans}.
#' @param xlev A list of levels for factors in the reference grid.
#' @param grid A data frame representing the reference grid for which EMMs are
#'   to be computed.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with the components required by \pkg{emmeans}:
#' \describe{
#'   \item{X}{Model matrix for the reference grid.}
#'   \item{bhat}{Estimated regression coefficients.}
#'   \item{V}{Variance-covariance matrix of the coefficients.}
#'   \item{nbasis}{A matrix indicating no penalized basis functions (empty in
#'     this case).}
#'   \item{dffun}{Function returning degrees of freedom for inference (fixed
#'     at \code{Inf}).}
#'   \item{dfargs}{Arguments for \code{dffun}.}
#'   \item{misc}{List of transformation and variance functions used in EMM
#'     computations.}
#' }
#'
#' @importFrom stats model.matrix vcov .getXlevels delete.response qnorm rnorm
#' @export
emm_basis.cpn <- function(object, trms, xlev, grid, ...) { # nolint

  # Construct the model matrix based on the reference grid and terms
  X <- model.matrix(trms, grid) # nolint

  # Extract regression coefficients
  bhat <- object$coefficients

  # Extract the full variance-covariance matrix of all estimated parameters
  V_full <- vcov(object) # nolint

  # Subset the variance-covariance matrix to include only regression parameters
  V <- V_full[names(bhat), names(bhat), drop = FALSE] # nolint

  # Ensure the column names of X match the coefficient names
  if (!all(names(bhat) %in% colnames(X))) {
    stop("Mismatch between model matrix columns and coefficient names")
  }

  # Reorder and subset X to match the order of bhat
  X <- X[, names(bhat), drop = FALSE] # nolint

  # Specify no penalized or nuisance basis functions (not used in this model)
  nbasis <- matrix(0, 0, length(bhat))
  colnames(nbasis) <- names(bhat)


  # Degrees of freedom used for inference
  dfargs <- list(df = Inf)

  # Simple function returning the residual df
  dffun <- function(k, dfargs) dfargs$df

  # Define inverse link function for transformation to response scale
  linkinv <- exp  # Used for log-link models

  # Define the variance function on the linear predictor scale
  varfun <- function(var) var^2  # Needed by emmeans for estimating SEs

  # Return the list required by emmeans to construct estimated marginal means
  list(
    X = X,
    bhat = bhat,
    V = V,
    nbasis = nbasis,
    dffun = dffun,
    dfargs = dfargs,
    misc = list(
      tran = "log",         # Indicates a log-link function
      inv.link = linkinv,   # How to transform back to response scale
      var.fun = varfun      # Variance function on the link scale
    )
  )
}
