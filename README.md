# smoothic <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/smoothic)](https://cran.r-project.org/package=smoothic)
![Download count](https://cranlogs.r-pkg.org/badges/grand-total/smoothic)
<!-- badges: end -->

For more information, check out the `smoothic`
[website](https://meadhbh-oneill.ie/smoothic/).

Implementation of the SIC epsilon-telescope method, either using single or multi-parameter regression. Includes classical regression with normally distributed errors and robust regression, where the errors are from the Laplace distribution. The "smooth generalized normal distribution" is used, where the estimation of an additional shape parameter allows the user to move smoothly between both types of regression. See O'Neill and Burke (2022) "Robust Distributional Regression with Automatic Variable Selection" for more details on [arXiv](https://arxiv.org/abs/2212.07317). This package also contains the data analyses from O'Neill and Burke (2023). "Variable selection using a smooth information criterion for distributional regression models" in [Statistics & Computing](https://link.springer.com/article/10.1007/s11222-023-10204-8). 

## Installation

### CRAN
You can install the released version of smoothic from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("smoothic")
```

### Github
Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("meadhbh-oneill/smoothic")
```

### smoothic()
The `smoothic()` function performs automatic variable selection with distributional regression.

``` r
library(smoothic)
fit <- smoothic(
  formula = y ~ .,
  data = dataset
)
```

A summary table that includes the estimated coefficients, estimated standard errors (SEE) and the value of the penalized likelihood function is returned with:

``` r
summary(fit)
```

Further information and examples of implementation (including plotting of the coefficient paths - `vignette("sgnd-boston")`) are available in the function documentation and vignettes.
