# Plot Diagnostics for a CPN Model

Generates diagnostic plots for a fitted Compound Poisson-Normal (CPN)
model object. Options include a residuals vs fitted values plot and a
Q-Q plot of the deviance residuals.

## Usage

``` r
# S3 method for class 'cpn'
plot(x, which = c("residuals", "qq"), ...)
```

## Arguments

- x:

  An object of class `"cpn"`, typically resulting from a CPN model
  fitting function.

- which:

  A character string specifying the type of plot to produce. Options are
  `"residuals"` (default) for a residuals vs fitted values plot, and
  `"qq"` for a Q-Q plot of deviance residuals.

- ...:

  Additional graphical parameters passed to the underlying plotting
  functions.

## Value

This function is called for its side effects and does not return a
value.

## Details

The residuals vs fitted plot helps assess non-linearity, unequal error
variances, and outliers. The Q-Q plot checks for normality of deviance
residuals.

## See also

[`residuals.cpn`](https://laszlopecze77.github.io/CPN/reference/residuals.cpn.md),
[`fitted.cpn`](https://laszlopecze77.github.io/CPN/reference/fitted.cpn.md)
