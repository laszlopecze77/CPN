# Extract coefficients from a Compound Poisson-Normal model

Returns the estimated regression coefficients from a fitted Compound
Poisson-Normal (CPN) model. Optionally includes the estimated
distribution parameters `mu` and `sigma`.

## Usage

``` r
# S3 method for class 'cpn'
coef(object, full = TRUE, ...)
```

## Arguments

- object:

  An object of class `cpn`.

- full:

  Logical. If `TRUE` (default), returns both the regression coefficients
  and the estimated distribution parameters `mu` and `sigma`. If
  `FALSE`, returns only the regression coefficients.

- ...:

  Additional arguments (currently ignored).

## Value

A named numeric vector of model coefficients. If `full = TRUE`, includes
`mu` and `sigma`.
