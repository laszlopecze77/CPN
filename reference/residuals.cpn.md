# Extract Residuals from a Compound Poisson-Normal (CPN) Model

This method extracts residuals from a fitted model object of class
`"cpn"`. Supports both raw and deviance residuals.

## Usage

``` r
# S3 method for class 'cpn'
residuals(object, type = c("deviance", "raw"), ...)
```

## Arguments

- object:

  An object of class `"cpn"`, typically produced by a CPN model fitting
  function.

- type:

  Character string indicating the type of residuals to return. Choices
  are `"deviance"` (default) or `"raw"`.

- ...:

  Additional arguments (currently unused).

## Value

A numeric vector of residuals. Returns deviance residuals if
`type = "deviance"`, or raw residuals (observed - fitted) if
`type = "raw"`.
