# Negative Log-Likelihood for Compound Poisson-Normal Regression

Computes the negative log-likelihood for a regression model where the
response is assumed to follow a Compound Poisson-Normal distribution.

## Usage

``` r
cpn_neg_log_likelihood(beta_mu_sigma, X, y, k_max = 10)
```

## Arguments

- beta_mu_sigma:

  Numeric vector. Contains the regression coefficients for the log-link
  Poisson mean (`lambda`), followed by the mean (`mu`) and standard
  deviation (`sigma`) of the Normal components.

- X:

  Numeric matrix. The design matrix for the regression model (rows =
  observations, columns = covariates).

- y:

  Numeric vector. The observed responses.

- k_max:

  Integer. Maximum number of Poisson events to consider for truncation
  in likelihood approximation (default is 10).

## Value

A numeric value giving the negative log-likelihood. Returns `Inf` if
invalid parameters are provided or numerical issues occur.

## Details

The function computes the likelihood for each observation as a weighted
sum over a finite number of Poisson-Normal mixtures, truncating at
`k_max`. If any resulting likelihood is zero, the function returns `Inf`
to penalize the parameter set.

## Examples

``` r
X <- matrix(c(1, 2, 3, 4), ncol = 2)
y <- c(0.5, 1.5)
beta_mu_sigma <- c(0.1, 0.2, 1, 0.5)
cpn_neg_log_likelihood(beta_mu_sigma, X, y, k_max = 10)
#> [1] 3.42991
```
