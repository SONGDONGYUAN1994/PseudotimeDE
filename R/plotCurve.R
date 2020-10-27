#' @title Plot the Fitted Model (Gene Trajectory)
#'
#' @description
#' Plot the fitted curve (gene trajectory) by a given gam model.
#' @param gene.vec A vector of genes. It should be a subset of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sce A SingleCellExperment object which contain the count data. Its row names should be genes and col names should be cells.
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
                      sce,
                      model.fit,
                      alpha = 0.2,
                      ncol = 2) {
  stopifnot(length(gene.vec) == length(model.fit))

  cell <- gene <- ori_pseudotime <- pseudotimes <- pseudotime <- counts <- NULL

  count_mat <- assays(sce)$counts

  count_mat <- cbind(t(count_mat), pseudotime = ori.tbl$pseudotime)

  count_mat <- count_mat %>% as_tibble() %>%
    tidyr::pivot_longer(cols = gene.vec, names_to = "gene", values_to = "counts") %>%
    dplyr::select(gene, pseudotime, counts)

  dat <- mapply(X = gene.vec, Y = model.fit, function(X, Y) {
    count_mat %>% dplyr::filter(gene ==  X) %>% dplyr::mutate(fitted = predict(Y, type = "response"))
  }, SIMPLIFY = FALSE)


  dat <- dplyr::bind_rows(dat)

  p <- dat %>% ggplot(aes(x = pseudotime, y = log10(counts+1))) + geom_point(alpha = alpha) +
    facet_wrap(~gene, ncol = ncol, scales = "free_y") +
    ylab("log10(count + 1)") +
    geom_line(aes(y = log10(fitted+1)), col = "blue", lty = "dashed", size = 1) +
    theme_bw()
  p
}

