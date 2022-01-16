# PseudotimeDE
`PseudotimeDE` is a robust method that accounts for the uncertainty in pseudotime inference and thus identifies DE genes along cell pseudotime with well-calibrated *p*-values. `PseudotimeDE` is flexible in allowing users to specify the pseudotime inference method and to choose the appropriate model for scRNA-seq data.

Latest News
------------
2021/12/02:
Update the vignettes

2021/11/12:
Replaced parallel.mclapply by BiocParallel.bplapply

2021/11/3:
Added QGAM(Smooth additive quantile regression model) as a model option.

2021/10/25:
Added Gaussian as a distribution option.

2021/10/7:
Added expression matrix and SeuratObject as input choices.

Introduction
------------
`PseudotimeDE` is developed to perfrom the differential expression (DE) test on genes along pseudotime (trajectory). Users can choose the pseudotime inference methods based on their preference. Basically, `PseudotimeDE` will use subsampling to capture the uncertainty of inferred pseudotime, and generate well-calibrated *p*-values.

Installation
------------

The package is not on Bioconductor or CRAN yet. For installation please use the following codes in `R`.

``` r
install.packages("devtools")
library(devtools)

devtools::install_github("SONGDONGYUAN1994/PseudotimeDE")
```
Please note that `PseudotimeDE` can be computationally intensive; we recommend users to allocate at least 10 cores unless they want to ignore the uncertainty of inferred pseudotime.

Quick start
-----------

For usage, please check the [vignettes](https://htmlpreview.github.io/?https://rpubs.com/dongyuansong/842884).
A separate [vignettes(QGAM)](http://rpubs.com/tianyang1998/855799) is created for PseudoimeDE when its model parameter is set as 'qgam'. 
If you meet problems, please contact <dongyuansong@ucla.edu>. 

Clarification
-----------
`PseudotimeDE` assumes that the latent pseudotime **does** exist and your pseudotime inference method can somehow capture it. If there is no latent pseudotime at all (e.g., the data is many iid random Poisson points), and you use the first PC of your data as the pseudotime, you will probabily fail to control the type 1 error (in a simpler word, you may still find some DE genes though the data is just random points). This phenomenon is related to the well-known **double dipping** problem. We greatly appreciate Dr. Daniela Witten, Dr. Anna Neufeld and Dr. Lucy Gao for raising this nice example.

This problem may not be very common on analysing real-world scRNA-seq data (since people usually believe that there is a latent structure), but it is interesting and important. We are now actively working on this problem!

Reference
-----------
Song, D., Li, J.J. PseudotimeDE: inference of differential gene expression along cell pseudotime with well-calibrated *p*-values from single-cell RNA sequencing data. *Genome Biol* **22**, 124 (2021). https://doi.org/10.1186/s13059-021-02341-y
