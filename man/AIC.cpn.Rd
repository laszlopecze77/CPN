% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AIC.R
\name{AIC.cpn}
\alias{AIC.cpn}
\title{Compute Akaike Information Criterion (AIC) for a CPN Model}
\usage{
\method{AIC}{cpn}(object, ..., k = 2)
}
\arguments{
\item{object}{An object of class \code{"cpn"}, typically produced by a CPN
model fitting function.}

\item{...}{Additional arguments (currently unused).}

\item{k}{Numeric penalty per parameter; the default is \code{2},
corresponding to the traditional AIC.}
}
\value{
A numeric value representing the AIC of the model.
}
\description{
Calculates the AIC for a fitted Compound Poisson-Normal (CPN) model.
}
\seealso{
\code{\link{logLik}}, \code{\link{BIC}}, \code{\link{AIC}}
}
