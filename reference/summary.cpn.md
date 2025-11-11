# Summarize a Compound Poisson-Normal (CPN) Model Fit

Produces a summary for a fitted `cpn` model object, including parameter
estimates, standard errors, z-values, p-values, deviance residual
summaries, and model diagnostics.

## Usage

``` r
# S3 method for class 'cpn'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"cpn"`, see
  [`cpn()`](https://laszlopecze77.github.io/CPN/reference/cpn.md).

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `summary.cpn`, which is a list containing:

- call:

  The matched call that generated the model.

- summary_table:

  A data frame with parameter estimates, standard errors, z-values, and
  p-values.

- deviance_summary:

  A five-number summary (Min, 1Q, Median, 3Q, Max) of deviance
  residuals.

- mu:

  Estimated value of the Normal component mean.

- sigma:

  Estimated value of the Normal component standard deviation.

- null_deviance:

  Null deviance of the model.

- residual_deviance:

  Residual deviance of the model.

- df_null:

  Degrees of freedom for the null model.

- df_residual:

  Degrees of freedom for the residual model.

- aic:

  Akaike Information Criterion for model comparison.

## See also

[`cpn()`](https://laszlopecze77.github.io/CPN/reference/cpn.md),
[`anova.cpn()`](https://laszlopecze77.github.io/CPN/reference/anova.cpn.md)

## Examples

``` r
set.seed(123)
data <- simulate_cpn_data()

# Sequential analysis of deviance
fit <- cpn(y ~ x1 + x2, data = data)
summary(fit)
#> Call:
#> cpn(formula = y ~ x1 + x2, data = data)
#> 
#> Deviance Residuals:
#>    Min     1Q Median     3Q    Max 
#> -2.871 -2.059 -1.269  2.181  3.745 
#> 
#> Coefficients:
#>             Estimate Std.Error z.value      Pr.z    
#> (Intercept)  0.63754   0.15190  4.1969 2.705e-05 ***
#> x1B         -0.60713   0.22934 -2.6473  0.008114 ** 
#> x2           0.53609   0.10791  4.9679 6.767e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Estimated mu parameter: 0.9600
#> Estimated sigma parameter: 1.7062
#> 
#> Null deviance: 478.20 on 99 degrees of freedom
#> Residual deviance: 449.74 on 95 degrees of freedom
#> AIC: 459.74
```
