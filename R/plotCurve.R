#' @title Plot the Fitted Model (Gene Trajectory)
#'
#' @description
#' Plot the fitted curve (gene trajectory) by a given gam model.
#' @param gene.vec A vector of genes. It should be a subset of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param mat The input expression data. It can be:
#' (1) A SingleCellExperment object which contain the expression data;
#' (2) An matrix;
#' (3) A Seurat object which contain the expression data.
#' Its row names should be genes and col names should be cells.
#' @param assay.use The \code{assay} used in SingleCellExperiment or \code{slot} used in Seurat. Default is \code{counts}.
#' @param model.fit A list of fitted gam models corresponding to genes in \code{gene.vec}.
#' @param alpha A numeric value of the opacity of points. Default is 0.2.
#' @param ncol A integer of facet column number. Default is 2.
#'
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assays colData
#' @importFrom mgcv gam
#' @importFrom tibble as_tibble
#'
#' @export plotCurve
#' @author Dongyuan Song

plotCurve <- function(gene.vec,
                      ori.tbl,
                      mat,
                      assay.use = "counts",
                      model.fit,
                      alpha = 0.2,
                      ncol = 2) {
  stopifnot(is.list(model.fit))
  stopifnot(length(gene.vec) == length(model.fit))


  log_counts <- cell <- gene <- ori_pseudotime <- pseudotimes <- pseudotime <- counts <- NULL

  # Subset mat
  mat <- mat[, ori.tbl$cell]

  if(class(mat)[1] == "SingleCellExperiment") {
    count_mat <- SummarizedExperiment::assay(mat, assay.use)
  }
  else if(class(mat)[1] == "SeuratObject") {
    count_mat <- Seurat::GetAssayData(object = mat, slot = assay.use)
  }
  else {
    count_mat <- mat
  }

  count_mat <- cbind(t(count_mat), pseudotime = ori.tbl$pseudotime)

  if(assay.use == "logcounts"){
    count_mat <- count_mat %>% as_tibble() %>%
      tidyr::pivot_longer(cols = gene.vec, names_to = "gene", values_to = "log_counts") %>%
      dplyr::select(gene, pseudotime, log_counts)
  }
  else{
    count_mat <- count_mat %>% as_tibble() %>%
      tidyr::pivot_longer(cols = gene.vec, names_to = "gene", values_to = "counts") %>%
      dplyr::select(gene, pseudotime, counts)
  }


  dat <- mapply(X = gene.vec, Y = model.fit, function(X, Y) {
    count_mat %>% dplyr::filter(gene ==  X) %>% dplyr::filter(!is.na(pseudotime)) %>% dplyr::mutate(fitted = predict(Y, type = "response"))
  }, SIMPLIFY = FALSE)


  dat <- dplyr::bind_rows(dat)

  if(assay.use == "logcounts"){
    p <- dat %>% ggplot(aes(x = pseudotime, y = log_counts)) + geom_point(alpha = alpha) +
      facet_wrap(~gene, ncol = ncol, scales = "free_y") +
      ylab("log10(count + 1)") +
      geom_line(aes(y = fitted), col = "blue", lty = "dashed", size = 1) +
      theme_bw() + theme(aspect.ratio=1)
  }
  else{
      p <- dat %>% ggplot(aes(x = pseudotime, y = log10(counts+1))) + geom_point(alpha = alpha) +
        facet_wrap(~gene, ncol = ncol, scales = "free_y") +
        ylab("log10(count + 1)") +
        geom_line(aes(y = log10(fitted+1)), col = "blue", lty = "dashed", size = 1) +
        theme_bw() + theme(aspect.ratio=1)
    }
  p
}

