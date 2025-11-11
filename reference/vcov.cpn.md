# Variance-Covariance Matrix for a CPN Model

Computes the variance-covariance matrix of parameter estimates from a
fitted Compound Poisson-Normal (CPN) regression model using the
numerical Hessian of the negative log-likelihood.

## Usage

``` r
# S3 method for class 'cpn'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `"cpn"` returned by the
  [`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md)
  function.

- ...:

  Additional arguments (currently unused).

## Value

A variance-covariance matrix of the model parameters. Rows and columns
are named according to the model parameters (regression coefficients,
`mu`, and `sigma`). If the Hessian is singular, contains `NA` values, or
is otherwise invalid, a matrix of `NA` values is returned with an
appropriate warning.

## Details

The Hessian matrix of the negative log-likelihood is computed using
numerical finite differences (via the
[`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html) function).
The variance-covariance matrix is then obtained by inverting this
Hessian. The result reflects local curvature and can be used to compute
standard errors and confidence intervals for the parameters.

## See also

[`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md),
[`summary.cpn`](https://laszlopecze77.github.io/CPN/reference/summary.cpn.md),
[`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)

## Examples

``` r
set.seed(123)
df <- data.frame(x = rnorm(100))
df$y <- sapply(exp(0.5 * df$x), function(lam) {
  k <- rpois(1, lam)
  if (k == 0) return(0)
  sum(rnorm(k, mean = 1, sd = 1))
})
fit <- cpn(y ~ x, data = df)
vcov(fit)
#>              (Intercept)             x           mu         sigma
#> (Intercept)  0.019229096 -0.0053012925 -0.007015153 -0.0042571726
#> x           -0.005301292  0.0169143475 -0.002521948  0.0008552595
#> mu          -0.007015153 -0.0025219477  0.016074768  0.0035950679
#> sigma       -0.004257173  0.0008552595  0.003595068  0.0145751698
```
