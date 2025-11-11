# Extract Log-Likelihood from a CPN Model

Returns the log-likelihood of a fitted Compound Poisson-Normal (CPN)
model.

## Usage

``` r
# S3 method for class 'cpn'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"cpn"`, typically resulting from a CPN model
  fitting function.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `"logLik"` representing the log-likelihood. The `df`
attribute contains the number of estimated parameters.

## See also

[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`BIC`](https://rdrr.io/r/stats/AIC.html),
[`logLik`](https://rdrr.io/r/stats/logLik.html)
