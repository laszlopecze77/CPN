# Basis for Estimated Marginal Means for Compound Poisson-Normal Models

Provides the required components for computing estimated marginal means
(EMMs) from a fitted Compound Poisson-Normal (CPN) model. This method is
used by the emmeans package to interface with models of class `"cpn"`.

## Usage

``` r
emm_basis.cpn(object, trms, xlev, grid, ...)
```

## Arguments

- object:

  A fitted model object of class `"cpn"` returned by a call to e.g.,
  `cpn5()`.

- trms:

  Terms object extracted from the model formula, typically provided by
  emmeans.

- xlev:

  A list of levels for factors in the reference grid.

- grid:

  A data frame representing the reference grid for which EMMs are to be
  computed.

- ...:

  Additional arguments (currently unused).

## Value

A list with the components required by emmeans:

- X:

  Model matrix for the reference grid.

- bhat:

  Estimated regression coefficients.

- V:

  Variance-covariance matrix of the coefficients.

- nbasis:

  A matrix indicating no penalized basis functions (empty in this case).

- dffun:

  Function returning degrees of freedom for inference (fixed at `Inf`).

- dfargs:

  Arguments for `dffun`.

- misc:

  List of transformation and variance functions used in EMM
  computations.
