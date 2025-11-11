# Update a CPN Model with a New Formula or Data

Updates a fitted Compound Poisson-Normal (CPN) model by changing the
model formula or data.

## Usage

``` r
# S3 method for class 'cpn'
update(object, new_formula, data = NULL, ...)
```

## Arguments

- object:

  An object of class `"cpn"`, typically produced by a CPN model fitting
  function.

- new_formula:

  A formula object specifying the new model formula.

- data:

  A data frame containing the new data. If `NULL`, the original data
  from the object is used.

- ...:

  Additional arguments passed to the underlying model fitting function.

## Value

A new `"cpn"` object fitted with the updated formula and/or data.

## See also

[`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md)
