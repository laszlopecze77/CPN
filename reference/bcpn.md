# Bayesian Compound Poisson-Normal (bCPN) Regression Model

Fits a Bayesian Compound Poisson-Normal regression model using JAGS.
Supports intercept-only or single binary factor models.

## Usage

``` r
bcpn(
  formula,
  data,
  n_iter = 20000,
  burn_in = 5000,
  n_chains = 1,
  k_max = 10,
  priors = list(beta_0 = "dnorm(0, 0.000001)", beta_1 = "dnorm(0, 0.00001)", mu =
    "dnorm(0, 0.000000001)", sigma = "dgamma(0.01, 0.01)")
)
```

## Arguments

- formula:

  A formula specifying the model. Must be either intercept-only
  (`y ~ 1`) or contain exactly one binary factor (`y ~ factor`).

- data:

  A data frame containing the variables in the model.

- n_iter:

  Integer. Number of MCMC iterations to run after burn-in. Default is
  20,000.

- burn_in:

  Integer. Number of MCMC burn-in iterations. Default is 5,000.

- n_chains:

  Integer. Number of MCMC chains. Default is 1.

- k_max:

  Numeric scalar or NULL. Maximum truncation value for the Poisson
  counts. Must be between 10 and 100. Default is 10.

- priors:

  A named list specifying prior distributions as JAGS model strings.
  Defaults are:

  - `beta_0`: prior for intercept (default: `"dnorm(0, 0.000001)"`)

  - `beta_1`: prior for factor levels (default: `"dnorm(0, 0.00001)"`)

  - `mu`: prior for normal mean (default: `"dnorm(0, 0.000000001)"`)

  - `sigma`: prior for normal sd (default: `"dgamma(0.01, 0.01)"`)

## Value

An `mcmc.list` object from the `coda` package containing posterior
samples for monitored parameters (including `mu`, `sigma`, `Expected`,
and `Diff` if applicable).

## Details

The model assumes the response follows a Compound Poisson-Normal
distribution where the count variable \\X\\ is Poisson-distributed with
a log-linear model for its rate parameter, and the observations \\y\\
are normal with mean \\\\\mu\\ X\\ and variance \\\\\sigma\\^2 / X\\
(with a large variance if \\X=0\\).

The parameter `k_max` truncates the Poisson rate at a maximum value to
avoid numerical issues in JAGS.

Only models with an intercept or a single binary factor are supported.

## Examples

``` r
testdata <- simulate_cpn_data()
fit <- bcpn(y ~ x1, data = testdata, n_iter = 20000)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 105
#>    Total graph size: 725
#> 
#> Initializing model
#> 
summary(fit)
#> 
#> Iterations = 6001:26000
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 20000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                Mean     SD Naive SE Time-series SE
#> Diff        -0.9023 0.3643 0.002576       0.004937
#> Expected[1]  2.0088 0.3657 0.002586       0.003215
#> Expected[2]  1.1065 0.2637 0.001864       0.003000
#> mu           1.1484 0.1979 0.001399       0.002527
#> sigma        1.9666 0.2127 0.001504       0.003129
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%     50%     75%   97.5%
#> Diff        -1.6545 -1.1354 -0.8893 -0.6575 -0.2135
#> Expected[1]  1.3355  1.7552  1.9926  2.2454  2.7692
#> Expected[2]  0.6599  0.9171  1.0860  1.2701  1.6836
#> mu           0.7785  1.0114  1.1413  1.2768  1.5571
#> sigma        1.5922  1.8175  1.9501  2.1020  2.4276
#> 
```
