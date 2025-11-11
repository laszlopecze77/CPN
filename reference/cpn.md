# Compound Poisson-Normal Regression

Fits a Compound Poisson-Normal (CPN) regression model to the response
variable using maximum likelihood estimation via the Nelder-Mead method
(Nelder and Mead 1965) .

## Usage

``` r
cpn(formula, data = NULL, mu_init = NULL, sigma_init = NULL, k_max = 10)
```

## Arguments

- formula:

  A formula specifying the model, e.g., `y ~ x1 + x2`.

- data:

  An optional data frame, list, or environment containing the variables
  in the model. If not found in `data`, the variables are taken from the
  environment from which `cpn` is called.

- mu_init:

  Optional initial value for the Normal mean parameter. If `NULL`, it is
  set to the sample mean of `y`.

- sigma_init:

  Optional initial value for the Normal standard deviation. If `NULL`,
  it is set to the sample standard deviation of `y`.

- k_max:

  Upper limit of summation used in approximating the Poisson
  convolution. Should be between 10 and 100. Default is 10.

## Value

An object of class `"cpn"` containing the following components:

- coefficients:

  Estimated regression coefficients.

- mu:

  Estimated mean of the Normal component.

- sigma:

  Estimated standard deviation of the Normal component.

- se:

  Standard errors of the estimated parameters.

- model:

  The model frame used.

- terms:

  The terms object from the model.

- formula:

  The model formula.

- data:

  The original data passed in.

- deviance_residuals:

  Vector of deviance residuals.

- fitted_values:

  Fitted values (Poisson mean times Normal mean).

- neg_log_likelihood:

  Negative log-likelihood at the optimum.

- null_deviance:

  Null model deviance.

- residual_deviance:

  Residual deviance of the fitted model.

- df_null:

  Degrees of freedom for the null model.

- df_residual:

  Degrees of freedom for the fitted model.

- k_max:

  Value used to truncate the Poisson convolution sum.

- call:

  The matched function call.

## Details

The function fits a regression model where the response is assumed to
follow a Compound Poisson-Normal distribution. The model estimates
regression coefficients for the Poisson intensity, and mean and standard
deviation of the Normal component.

The optimization is performed via maximum likelihood using the
Nelder-Mead method. Standard errors are computed via the observed
information matrix using numerical Hessian. If the Hessian is singular
or not positive definite, a warning is issued and standard errors are
returned as `NA`.

## References

Nelder JA, Mead R (1965). “A simplex method for function minimization.”
*The computer journal*, **7**(4), 308–313.
[doi:10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308) .

## Examples

``` r
set.seed(123)
data <- simulate_cpn_data()

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
