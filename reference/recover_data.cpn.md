# Recover Data for `cpn` Model Objects

This method is used by the emmeans package to recover the data and terms
from a fitted Compound Poisson-Normal (CPN) model object. It is
typically used internally when computing estimated marginal means.

## Usage

``` r
# S3 method for class 'cpn'
recover_data(object, ...)
```

## Arguments

- object:

  A fitted model object of class `cpn`.

- ...:

  Additional arguments (not used).

## Value

An object containing the recovered data and terms, suitable for use by
emmeans.
