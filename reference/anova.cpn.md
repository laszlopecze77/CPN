# Analysis of Deviance for CPN Models

Computes an analysis of deviance table for objects of class `cpn`,
either:

1.  sequentially by model terms (Type I ANOVA) when one model is
    supplied, or

2.  by comparing two or more nested models.

## Usage

``` r
# S3 method for class 'cpn'
anova(object, ...)
```

## Arguments

- object:

  An object of class `cpn`, typically the result of a call to
  [`cpn()`](https://laszlopecze77.github.io/CPN/reference/cpn.md).

- ...:

  Optional additional objects of class `cpn` for model comparison.

## Value

An object of class `"anova"` inheriting from `"data.frame"`, containing:

- If one model is provided: sequential deviance table by added terms.

- If multiple models are provided: deviance comparison between models.

## Details

When a single model is supplied, this function refits the model with
terms added sequentially and reports the change in residual deviance.

When multiple models are supplied, they are assumed to be nested and
ordered by increasing complexity (based on residual degrees of freedom).
The function then compares their deviances using chi-squared tests.

## Examples

``` r
# Simulated data
set.seed(123)
data <- simulate_cpn_data()

# Sequential analysis of deviance
fit <- cpn(y ~ x1 + x2, data = data)
anova(fit)
#>       Term     Df Deviance Resid. Df Resid. Dev   Pr(>Chi) Signif
#>  Residuals                        97      478.2                  
#>         x1      1   6.7993        96      471.4  0.0091193     **
#>         x2      1   21.657        95     449.74 3.2607e-06    ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Model comparison
fit1 <- cpn(y ~ x1, data = data)
fit2 <- cpn(y ~ x1 + x2, data = data)
anova(fit1, fit2)
#>    Model Resid. Df Resid. Dev     Df Deviance   Pr(>Chi) Signif
#>  Model 1        95     449.74                                  
#>  Model 2        96      471.4      1   21.657 3.2607e-06    ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
