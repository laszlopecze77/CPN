# Predict Method for CPN Model Objects

Computes predictions from a fitted Compound Poisson-Normal (CPN)
regression model. Supports predictions on the link, rate, or response
scale, with optional confidence intervals.

## Usage

``` r
# S3 method for class 'cpn'
predict(
  object,
  newdata = NULL,
  type = c("link", "rate", "response"),
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  An object of class `cpn`, typically the result of a call to a function
  fitting a Compound Poisson-Normal regression model.

- newdata:

  An optional data frame in which to look for variables with which to
  predict. If omitted, the original model data is used.

- type:

  Type of prediction: `"link"` returns the linear predictor \\\eta =
  X\beta\\; `"rate"` returns \\\exp(\eta)\\; `"response"` returns the
  mean response \\E\[Y\] = \mu \cdot \exp(\eta)\\.

- interval:

  Type of interval calculation. Either `"none"` (default) or
  `"confidence"` for confidence intervals around the predicted values.

- level:

  Confidence level for the interval. Defaults to 0.95.

- ...:

  Further arguments passed to or from other methods (not currently
  used).

## Value

If `interval = "none"`, returns a numeric vector of predicted values on
the specified scale. If `interval = "confidence"`, returns a data frame
with columns:

- `fit`:

  Predicted value

- `lwr`:

  Lower bound of the confidence interval

- `upr`:

  Upper bound of the confidence interval

## Details

For predictions on the response scale with confidence intervals, the
standard errors of both the linear predictor and the estimated `mu`
parameter are combined using the delta method.

Factor levels in `newdata` are aligned to match those used in the
original model fit.

## See also

[`cpn`](https://laszlopecze77.github.io/CPN/reference/cpn.md),
[`vcov`](https://rdrr.io/r/stats/vcov.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`predict`](https://rdrr.io/r/stats/predict.html)

## Examples

``` r
set.seed(123)
data <- simulate_cpn_data()

fit <- cpn(y ~ x1 + x2, data = data)
predict(fit, type = "response", interval = "confidence")
#>           fit        lwr       upr
#> 1   2.0802264 1.14759018 3.0128627
#> 2   1.7884887 0.98146594 2.5955114
#> 3   1.7748077 0.97332482 2.5762907
#> 4   2.0611074 0.88766658 3.2345483
#> 5   1.6090445 0.87228376 2.3458053
#> 6   2.2311441 0.92886068 3.5334275
#> 7   0.4313962 0.15940715 0.7033853
#> 8   1.3538511 0.65162882 2.0560733
#> 9   1.9407450 1.06996226 2.8115278
#> 10  2.0389587 1.12496774 2.9529497
#> 11  1.2129653 0.59009048 1.8358402
#> 12  0.7559761 0.35337806 1.1585741
#> 13  0.8277172 0.39413091 1.2613035
#> 14  1.0519275 0.50624821 1.5976067
#> 15  0.5570876 0.23539798 0.8787771
#> 16  2.1369808 1.17822896 3.0957327
#> 17  1.2583836 0.61051292 1.9062543
#> 18  1.8684137 1.02840620 2.7084212
#> 19  2.9775332 1.57060806 4.3844582
#> 20  5.4505108 2.21924635 8.6817753
#> 21  1.3957565 0.73621730 2.0552958
#> 22  0.2869689 0.07575121 0.4981865
#> 23  3.1137985 1.62418864 4.6034083
#> 24  1.2416952 0.63439214 1.8489982
#> 25  1.2558825 0.64387032 1.8678946
#> 26  3.1470817 1.63688829 4.6572752
#> 27  0.8494906 0.40625913 1.2927221
#> 28  0.5143398 0.20950267 0.8191769
#> 29  2.0014462 1.10415223 2.8987402
#> 30  0.9185928 0.44395795 1.3932276
#> 31  1.8216902 1.00109458 2.6422859
#> 32  1.2166389 0.59176324 1.8415147
#> 33  1.4887943 0.79634817 2.1812405
#> 34  1.3979286 0.66979960 2.1260576
#> 35  0.8792773 0.42265961 1.3358950
#> 36  2.1695948 1.19558811 3.1436014
#> 37  3.2696457 1.68237893 4.8569124
#> 38  2.2932547 1.25977745 3.3267320
#> 39  1.5249250 0.81938413 2.2304658
#> 40  1.8320090 0.82334289 2.8406751
#> 41  3.0934422 1.61634729 4.5705370
#> 42  1.3278188 0.64065793 2.0149798
#> 43  1.1247137 0.54879019 1.7006371
#> 44  1.2970066 0.67123719 1.9227760
#> 45  3.7663637 1.84745526 5.6852721
#> 46  1.3163728 0.68406598 1.9486796
#> 47  5.8666695 2.27493524 9.4584038
#> 48  2.2505329 0.93322569 3.5678401
#> 49  1.6005022 0.86696001 2.3340444
#> 50  1.0475124 0.50324662 1.5917783
#> 51  0.6761791 0.30677299 1.0455852
#> 52  2.0842061 1.14975646 3.0186558
#> 53  1.5910991 0.86108694 2.3211112
#> 54  1.5073599 0.80820802 2.2065118
#> 55  1.0903723 0.53236148 1.6483831
#> 56  0.9659987 0.46909608 1.4629013
#> 57  0.6497062 0.29106153 1.0083508
#> 58  0.7426754 0.29727708 1.1880736
#> 59  0.8071140 0.38254885 1.2316792
#> 60  2.9723167 1.56850535 4.3761280
#> 61  1.3340715 0.69575521 1.9723878
#> 62  1.3709052 0.65871940 2.0830910
#> 63  0.4157013 0.15001446 0.6813881
#> 64  1.7627733 0.96613798 2.5594085
#> 65  2.3991748 1.31273035 3.4856193
#> 66  1.1629876 0.56696554 1.7590097
#> 67  1.9219241 1.05923431 2.7846139
#> 68  1.2881371 0.66534883 1.9109253
#> 69  1.1516027 0.57382304 1.7293824
#> 70  1.0488004 0.50412229 1.5934785
#> 71  1.0540245 0.51415463 1.5938943
#> 72  1.0927973 0.53400691 1.6515877
#> 73  1.3961110 0.73644853 2.0557735
#> 74  1.5831010 0.85608099 2.3101210
#> 75  4.8800474 2.12156838 7.6385264
#> 76  0.6977048 0.31946307 1.0759466
#> 77  1.1226985 0.54782200 1.6975750
#> 78  1.8935793 1.04296462 2.7441940
#> 79  0.5909063 0.25582175 0.9259909
#> 80  0.9524844 0.46199122 1.4429777
#> 81  2.1467582 0.90908086 3.3844355
#> 82  1.2606079 0.61149869 1.9097172
#> 83  1.8566602 1.02157024 2.6917502
#> 84  0.7890299 0.37230057 1.2057592
#> 85  0.3291690 0.09935721 0.5589809
#> 86  1.8149310 0.81811541 2.8117465
#> 87  0.8299720 0.35561386 1.3043301
#> 88  2.7002783 1.45336967 3.9471868
#> 89  2.7538558 1.02496824 4.4827433
#> 90  0.8374569 0.36065538 1.3142585
#> 91  1.4416199 0.68731561 2.1959242
#> 92  0.8598341 0.41197957 1.3076886
#> 93  0.7818129 0.32331323 1.2403126
#> 94  0.4393515 0.16418328 0.7145197
#> 95  0.4193602 0.15220024 0.6865202
#> 96  1.3662363 0.71690917 2.0155634
#> 97  0.8294758 0.35527983 1.3036718
#> 98  1.4309422 0.68308000 2.1788045
#> 99  5.5986581 2.24050341 8.9568129
#> 100 0.9109270 0.41036220 1.4114917
```
