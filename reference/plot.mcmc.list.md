# Plot Posterior Distributions of Expected Means from bCPN Model

Plots posterior density estimates of the `Expected` parameter(s) from an
MCMC output of the `bcpn` function.

## Usage

``` r
# S3 method for class 'mcmc.list'
plot(
  x,
  param_names = NULL,
  col = c("skyblue", "salmon"),
  lwd = 2,
  add_medians = TRUE,
  ...
)
```

## Arguments

- x:

  An `mcmc.list` object containing samples from `bcpn`.

- param_names:

  Character vector of parameter names to plot. Defaults to all
  `Expected` parameters found.

- col:

  Colors for the density lines. Default is `c("skyblue", "salmon")`.

- lwd:

  Line width for density lines. Default is 2.

- add_medians:

  Logical indicating whether to add median vertical lines. Default TRUE.

- ...:

  Additional graphical parameters passed to `plot`.
