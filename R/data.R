#' A SingelCellExperiment dataset containing 4,016 genes and 390 cells
#' @format A dataset with 4,016 rows and 390 cols
"LPS_sce"

#' A tibble containing the inferred pseudotime of \code{LPS_sce} by \code{slingshot} package
#' @format A tibble with two cols
#' \describe{
#' \item{cell}{Cell names}
#' \item{pseudotime}{Inferred pseudotime of each cell}
#' }
"LPS_ori_tbl"

#' A list of tibbles containing the inferred pseudotime of each subsample
#' @format A list with 1,000 tibbles
"LPS_sub_tbl"
