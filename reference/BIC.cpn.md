# Bayesian Information Criterion for Compound Poisson-Normal Models

Computes the Bayesian Information Criterion (BIC) for a fitted Compound
Poisson-Normal (CPN) regression model.

## Usage

``` r
BIC.cpn(object, ...)
```

## Arguments

- object:

  An object of class `"cpn"`, typically the result of a call to
  [`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md).

- ...:

  Additional arguments (currently unused).

## Value

A numeric value representing the BIC of the fitted model.

## Details

The BIC is computed as: \$\$-2 \cdot \log L + k \cdot \log(n)\$\$ where
\\L\\ is the likelihood of the fitted model, \\k\\ is the number of
estimated parameters (including regression coefficients, \\\mu\\, and
\\\sigma\\), and \\n\\ is the number of observations.

## See also

[`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md),
[`AIC`](https://rdrr.io/r/stats/AIC.html),
[`BIC`](https://rdrr.io/r/stats/AIC.html)
