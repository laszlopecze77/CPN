#' Basis for Estimated Marginal Means for `cpn` Model
#'
#' This function provides the required components for computing estimated marginal means
#' from a fitted Compound Poisson-Normal (CPN) model using the \pkg{emmeans} package.
#' It extracts the model matrix, fixed-effect estimates, and covariance matrix, and
#' provides information about degrees of freedom and transformation on the response scale.
#'
#' @param object A fitted model object of class \code{cpn}.
#' @param trms Terms object, usually derived from \code{delete.response(terms(object))}.
#' @param xlev A list of levels for factors, typically from the data used in fitting.
#' @param grid A reference grid (data frame) over which to compute estimated means.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with components used by \pkg{emmeans}:
#' \describe{
#'   \item{\code{X}}{Model matrix for the reference grid.}
#'   \item{\code{bhat}}{Estimated regression coefficients.}
#'   \item{\code{V}}{Covariance matrix for the fixed effects.}
#'   \item{\code{nbasis}}{Matrix for any non-estimable basis (empty here).}
#'   \item{\code{dffun}}{Function to return degrees of freedom.}
#'   \item{\code{dfargs}}{Arguments to \code{dffun}.}
#'   \item{\code{misc}}{Miscellaneous information including link function details.}
#' }
#'
#' @importFrom emmeans emm_basis
#' @export
emm_basis.cpn <- function(object, trms, xlev, grid, ...) {
  # Construct model matrix from reference grid
  X <- model.matrix(trms, grid)

  # Extract regression coefficients only (exclude mu and sigma)
  bhat <- object$coefficients  # Named vector

  # Extract full vcov and reduce to regression part
  V_full <- vcov(object)
  V <- V_full[names(bhat), names(bhat), drop = FALSE]

  # Match X columns with coefficient names
  if (!all(names(bhat) %in% colnames(X))) {
    stop("Mismatch between model matrix columns and coefficient names")
  }
  X <- X[, names(bhat), drop = FALSE]

  # No penalized/nuisance terms
  nbasis <- matrix(0, 0, length(bhat))
  colnames(nbasis) <- names(bhat)

  # Degrees of freedom
  dfargs <- list(df = if (!is.null(object$df_residual)) object$df_residual else Inf)
  dffun <- function(k, dfargs) dfargs$df

  # Specify transformation for emmeans on response scale
  linkinv <- exp  # Inverse of log link
  varfun <- function(var) var^2  # Variance function on log scale (used by emmeans)

  list(
    X = X,
    bhat = bhat,
    V = V,
    nbasis = nbasis,
    dffun = dffun,
    dfargs = dfargs,
    misc = list(
      tran = "log",          # Tell emmeans that this is a log-link model
      inv.link = linkinv,   # How to go back to response scale
      var.fun = varfun      # Variance function
    )
  )
}
