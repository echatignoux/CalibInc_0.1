CalibInc
========

Code and material for the article *"For a sound use of healthcare in epidemiology: evaluation of a calibration model for count data with application to prediction of cancer incidence in areas without cancer registry."*

R-package
---------

The material and R code used in the paper are made available as an R-package named `CalibInc`. Use [devtools](https://github.com/hadley/devtools) to install the package from Github:

``` r
require(devtools)
install_github("echatignoux/CalibInc")
```

After installation, the package can be loaded into R.

``` r
library(CalibInc)
```

The package needs R version 3.1.2 or higher.

Analysis from the application section of the paper
--------------------------------------------------

The use of `CalibInc` package is illustrated with the code used in the
application section of the paper. It illustrates how to predict cancer
incidence for LOP cancer in men using registry and hospitalization
data. Code and results are available at
<https://rawgit.com/echatignoux/CalibInc/master/Application/Application.html>.
