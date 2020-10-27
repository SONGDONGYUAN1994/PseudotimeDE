#' @title Perform Differential Expression Test on a List of Genes
#'
#' @description
#' Test if each gene in a list of genes is differentially expressed along pseudotime.
#' A wrapper of \code{PseudotimeDE::pseudotimeDE}.
#'
#' @param gene.vec A vector of genes. It should be a subset of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sub.tbl A list of tibbles or dataframes where each is the fit of a subsample. Each element is the same format as ori.tbl.
#' @param sce A SingleCellExperment object which contain the count data. Its row names should be genes and col names should be cells.
#' @param model A string of the model name. One of \code{nb}, \code{zinb} and \code{auto}.
#' @param k A integer of the basis dimension. Default is 6. The reults are usually robust to different k; we recommend to use k from 5 to 10.
#' @param knots A numeric vector of the location of knots.
#' @param fix.weight A logic variable indicating if the ZINB-GAM will use the zero weights from the original model.
#' @param aicdiff A numeric variable of the threshold of model selection. Only works when \code{model = `auto`}.
#' @param seed A numeric variable of the random seed. It mainly affects the fitting of null distribution.
#' @param mc.cores Number of cores for computing.
#' @param SIMPLIFY A logic variable whether to return a tibble (TRUE) or a list of lists (FALSE). Default is TRUE.
#' @return A tibble of summary results of genes
#'
#' @examples
#' data("LPS_sce")
#' data("LPS_ori_tbl")
#' data("LPS_sub_tbl")
#' res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10"),
#' ori.tbl = LPS_ori_tbl,
#' sub.tbl = LPS_sub_tbl,
#' sce = LPS_sce,
#' model = "nb")
#'
#' @export runPseudotimeDE
#' @author Dongyuan Song

runPseudotimeDE <- function(gene.vec,
                            ori.tbl,
                            sub.tbl,
                            sce,
                            model = c("nb", "zinb", "auto"),
                            k = 6,
                            knots = c(0:5/5),
                            fix.weight = TRUE,
                            aicdiff = 10,
                            seed = 123,
                            mc.cores = 2,
                            SIMPLIFY = TRUE) {
  set.seed(seed)
  res <- parallel::mclapply(gene.vec, PseudotimeDE::pseudotimeDE,
                            ori.tbl = ori.tbl,
                            sub.tbl = sub.tbl,
                            sce = sce,
                            model = model,
                            k = k,
                            knots = knots,
                            fix.weight = fix.weight,
                            aicdiff = aicdiff,
                            mc.cores = mc.cores)
  # Avoid package check error
  expv.quantile <- gam.fit <- NULL

  if(SIMPLIFY) {
    res <- simplify2array(res)
    res <- t(res)
    rownames(res) <- gene.vec
    res <- tibble::as_tibble(res, rownames = "gene")
    res <- tidyr::unnest(res, cols = ! (gam.fit | expv.quantile))
    }

  res
}
