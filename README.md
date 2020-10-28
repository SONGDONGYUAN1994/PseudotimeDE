# PseudotimeDE
`PseudotimeDE` is a robust method that accounts for the uncertainty in pseudotime inference and thus identifies DE genes along cell pseudotime with well-calibrated *p*-values. `PseudotimeDE` is flexible in allowing users to specify the pseudotime inference method and to choose the appropriate model for scRNA-seq data.

Introduction
------------
`PseudotimeDE` is developed to perfrom the differential expression (DE) test on genes along pseudotime (trajectory). Users can choose the pseudotime inference methods based on their preference. Basically, `PseudotimeDE` will use subsampling to capture the uncertainty of inferred pseudotime, and generate well-calibrated *p*-values.

Installation
------------

The package is not on Bioconductor or CRAN yet. For installation please use the following codes in `R`.

``` r
install.packages("devtools")
library(devtools)

install_github("SONGDONGYUAN1994/PseudotimeDE")
```
Please notice that `PseudotimeDE` can be computationally intensive; we recommend users to allocate at least 10 cores unless they want to ignore the uncertainty of inferred pseudotime.

Quick start
-----------

For usage, please check the [Vignette](https://github.com/SONGDONGYUAN1994/PseudotimeDE/blob/master/vignettes/quickstart.Rmd).
If you meet problems, please contact <dongyuansong@ucla.edu>. 

