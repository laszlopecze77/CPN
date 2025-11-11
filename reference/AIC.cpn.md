# Compute Akaike Information Criterion (AIC) for a CPN Model

Calculates the AIC for a fitted Compound Poisson-Normal (CPN) model.

## Usage

``` r
# S3 method for class 'cpn'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  An object of class `"cpn"`, typically produced by a CPN model fitting
  function.

- ...:

  Additional arguments (currently unused).

- k:

  Numeric penalty per parameter; the default is `2`, corresponding to
  the traditional AIC.

## Value

A numeric value representing the AIC of the model.

## See also

[`logLik`](https://rdrr.io/r/stats/logLik.html),
[`BIC`](https://rdrr.io/r/stats/AIC.html),
[`AIC`](https://rdrr.io/r/stats/AIC.html)
