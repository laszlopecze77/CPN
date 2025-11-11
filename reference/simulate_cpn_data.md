# Simulate Data from a Compound Poisson-Normal Model

This function generates synthetic data based on a Compound
Poisson-Normal (CPN) model. The number of events for each observation is
drawn from a Poisson distribution, and the outcome is the sum of
normally distributed values for each event.

## Usage

``` r
simulate_cpn_data(n = 100, beta = c(0.5, -0.3, 0.7), mu = 1, sigma = 2)
```

## Arguments

- n:

  Integer. Number of observations to simulate. Default is 100.

- beta:

  Numeric vector. Coefficients for the linear predictor, including
  intercept. Default is `c(0.5, -0.3, 0.7)`. Must match the number of
  columns in the model matrix: intercept, x1B, x2.

- mu:

  Numeric. Mean of the Normal distribution for each event. Default is 1.

- sigma:

  Numeric. Standard deviation of the Normal distribution for each event.
  Default is 2.

## Value

A `data.frame` with `n` rows and 3 columns:

- y:

  Numeric response variable generated from the CPN model.

- x1:

  Categorical predictor with levels "A" and "B".

- x2:

  Continuous predictor drawn from a standard normal distribution.

## Details

The predictors include a binary categorical variable (`x1`) and a
continuous variable (`x2`). A linear model with coefficients `beta` is
used to model the log of the Poisson rate.

## Examples

``` r
set.seed(42)
simulated_data <- simulate_cpn_data(n = 200)
head(simulated_data)
#>            y x1          x2
#> 1  0.0000000  A -0.71040656
#> 2 -2.4874819  A  0.25688371
#> 3  0.3583388  A -0.24669188
#> 4  0.0000000  B -0.34754260
#> 5  0.0000000  A -0.95161857
#> 6 16.0883605  B -0.04502772
```
