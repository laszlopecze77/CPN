% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emm_basis.R
\name{emm_basis.cpn}
\alias{emm_basis.cpn}
\title{Basis for Estimated Marginal Means for Compound Poisson-Normal Models}
\usage{
emm_basis.cpn(object, trms, xlev, grid, ...)
}
\arguments{
\item{object}{A fitted model object of class \code{"cpn"} returned by a call to
e.g., \code{cpn5()}.}

\item{trms}{Terms object extracted from the model formula, typically
provided by \pkg{emmeans}.}

\item{xlev}{A list of levels for factors in the reference grid.}

\item{grid}{A data frame representing the reference grid for which EMMs are
to be computed.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A list with the components required by \pkg{emmeans}:
\describe{
\item{X}{Model matrix for the reference grid.}
\item{bhat}{Estimated regression coefficients.}
\item{V}{Variance-covariance matrix of the coefficients.}
\item{nbasis}{A matrix indicating no penalized basis functions (empty in
this case).}
\item{dffun}{Function returning degrees of freedom for inference (fixed
at \code{Inf}).}
\item{dfargs}{Arguments for \code{dffun}.}
\item{misc}{List of transformation and variance functions used in EMM
computations.}
}
}
\description{
Provides the required components for computing estimated marginal means
(EMMs) from a fitted Compound Poisson-Normal (CPN) model. This method is
used by the \pkg{emmeans} package to interface with models of class \code{"cpn"}.
}
