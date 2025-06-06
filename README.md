# CPN

## Installation

Install the current development version from GitHub with:

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("laszlopecze77/CPN")
```

## Development

Development dependencies:

- R version >= 4.5.0
- For windows machines, rtools45
  (<https://cran.rstudio.com/bin/windows/Rtools/rtools45/rtools.html>)
- package`roxygen2` version 7.3.2 for the documentation.
- package `Rdpack` for the citation within function documentation.
- package `devtools` to ease the development.
- package `rmarkdown` for the vignettes.
- package `sessioninfo` version 1.2.2 if running under windows.
  See current issue <https://github.com/r-lib/sessioninfo/issues/112>.
