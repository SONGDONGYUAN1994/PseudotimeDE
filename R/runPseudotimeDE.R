#' @title Perform Differential Expression Test on a List of Genes
#'
#' @description
#' Test if each gene in a list of genes is differentially expressed along pseudotime.
#' A wrapper of \code{PseudotimeDE::pseudotimeDE}.
#'
#' @param gene.vec A vector of genes. It should be a subset of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sub.tbl A list of tibbles or dataframes where each is the fit of a subsample. Each element is the same format as ori.tbl.
#' @param mat The input expression data. It can be:
#' (1) A SingleCellExperment object which contain the expression data;
#' (2) An matrix;
#' (3) A Seurat object which contain the expression data.
#' Its row names should be genes and col names should be cells.
#' @param assay.use The \code{assay} used in SingleCellExperiment or \code{slot} used in Seurat. Default is \code{counts}.
#' @param model A string of the model name. One of \code{nb}, \code{zinb}, \code{gaussian}, \code{auto} and \code{qgam}.
#' @param k A integer of the basis dimension. Default is 6. The reults are usually robust to different k; we recommend to use k from 5 to 10.
#' @param knots A numeric vector of the location of knots.
#' @param fix.weight A logic variable indicating if the ZINB-GAM will use the zero weights from the original model.
#' @param aicdiff A numeric variable of the threshold of model selection. Only works when \code{model = `auto`}.
#' @param seed A numeric variable of the random seed. It mainly affects the fitting of null distribution.
#' @param quant The quantile of interest for quantile regression (qgam), range from 0 to 1, default as 0.5.
#' @param usebam A logical variable. If use \code{mgcv::bam}, which may be faster with large sample size (e.g., > 10'000 cells).
#' @param mc.cores Number of cores for computing.
#' @param mc.preschedule See \code{mclapply}. Default is TRUE.
#' @param SIMPLIFY A logic variable whether to return a tibble (TRUE) or a list of lists (FALSE). Default is TRUE.
#' @return A tibble of summary results of genes
#' @importFrom BiocParallel bpparam bplapply
#'
#' @examples
#' data("LPS_sce")
#' data("LPS_ori_tbl")
#' data("LPS_sub_tbl")
#' res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10"), ori.tbl = LPS_ori_tbl, sub.tbl = LPS_sub_tbl[1:10], mat = LPS_sce, model = "nb")
#'
#' @export runPseudotimeDE
#' @author Dongyuan Song

runPseudotimeDE <- function(gene.vec,
                            ori.tbl,
                            sub.tbl,
                            mat,
                            assay.use = "counts",
                            model = c("nb", "zinb", "gaussian","auto", "qgam"),
                            k = 6,
                            knots = c(0:5/5),
                            fix.weight = TRUE,
                            aicdiff = 10,
                            seed = 123,
                            quant = 0.5,
                            usebam = FALSE,
                            mc.cores = 2,
                            mc.preschedule = TRUE,
                            SIMPLIFY = TRUE) {
  set.seed(seed)

  # Avoid package check error
  expv.quantile <- gam.fit <- NULL

  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- mc.cores

  res <- BiocParallel::bplapply(gene.vec, function(x, ...) {

    cur_res <- tryCatch(expr = PseudotimeDE::pseudotimeDE(gene = x, ...), error = function(e) {
      return(list(fix.pv = NA,
                  emp.pv = NA,
                  para.pv = NA,
                  ad.pv = NA,
                  rank = NA,
                  gam.fit = NA,
                  zinf = NA,
                  aic = NA,
                  expv.quantile = NA,
                  expv.mean = NA,
                  expv.zero = NA
      ))
    })
    cur_res
  },
                            ori.tbl = ori.tbl,
                            sub.tbl = sub.tbl,
                            mat = mat,
                            assay.use = assay.use,
                            model = model,
                            k = k,
                            knots = knots,
                            fix.weight = fix.weight,
                            aicdiff = aicdiff,
                            quant = quant,
                            usebam = usebam,
                            BPPARAM = BPPARAM)


  if(SIMPLIFY) {
    res <- simplify2array(res)
    res <- t(res)
    rownames(res) <- gene.vec
    res <- tibble::as_tibble(res, rownames = "gene")
    res <- tidyr::unnest(res, cols = ! (gam.fit | expv.quantile))
    }

  res
}
