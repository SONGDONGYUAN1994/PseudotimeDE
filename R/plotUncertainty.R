#' @title Plot the Pseudotime of Subsamples
#'
#' @description
#' Plot the fitted pseudotime of each cell from all subsamples.
#'
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sub.tbl A list of tibbles or dataframes where each is the fit of a subsample. Each element is the same format as \code{ori.tbl}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export plotUncertainty
#' @author Dongyuan Song, Huy Nyugen

plotUncertainty <- function(ori.tbl,
                            sub.tbl) {

  n_subample <- length(sub.tbl)
  cell <- gene <- ori_pseudotime <- pseudotimes <- pseudotime <- counts <- density <- NULL

  Cells_true_time <- as.data.frame(ori.tbl)
  colnames(Cells_true_time) <- c("cell", "ori_pseudotime")

  #sub.list <- lapply(sub.tbl, function(x) {x$pseudotime})

  Merge_pseudotimes <- suppressWarnings(Reduce(function(x,y) merge(x = x, y = y, by = "cell", all = TRUE), sub.tbl))

  for(i in 2:length(sub.tbl)+1){
    colnames(Merge_pseudotimes)[i] <- paste("Pseudotime", i-1, sep = " ")
  }

  Truetimes_pseudotimes <- merge(Cells_true_time, Merge_pseudotimes, by = "cell")
  Descending_truetimes_pseudotimes <- dplyr::arrange(Truetimes_pseudotimes, dplyr::desc(ori_pseudotime))

  pseudotimes_only <- Descending_truetimes_pseudotimes %>% dplyr::select(c(-cell, -ori_pseudotime))
  Descending_pseudotimes_transpose <- purrr::transpose(pseudotimes_only)
  Descending_unlist_pseudotimes <- data.frame(pseudotimes = unlist(Descending_pseudotimes_transpose, use.names = FALSE))

  cell_truetime_only <- Descending_truetimes_pseudotimes %>% dplyr::select(cell, ori_pseudotime)
  multiple_cell_truetime <- cell_truetime_only[rep(seq_len(nrow(cell_truetime_only)), each = n_subample), ]

  na_plotdata <- cbind(Descending_unlist_pseudotimes,multiple_cell_truetime)

  plotdata <- na_plotdata %>% stats::na.omit(pseudotimes)
  plotdata$cell <- factor(plotdata$cell, levels = unique(plotdata$cell))

  p <- ggplot(plotdata, aes(pseudotimes, cell)) +
    stat_density(aes(fill = after_stat(density)), geom = "raster", position = "identity") +
    scale_fill_gradient(low = "white", high = "black") +
    labs(x="Pseudotimes of subsamples", y = "Cells", fill = "Density") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y=element_blank(),
          aspect.ratio = 1, legend.position = "right")

  p
}



