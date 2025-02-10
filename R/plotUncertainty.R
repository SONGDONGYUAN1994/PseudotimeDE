#' @title Plot the Pseudotime of Subsamples
#'
#' @description
#' Plot the fitted pseudotime of each cell from all subsamples.
#'
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sub.tbl A list of tibbles or dataframes where each is the fit of a subsample. Each element is the same format as \code{ori.tbl}.
#'
#' @return A ggplot object visualizing pseudotime distributions across subsamples.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr arrange desc select slice mutate
#' @importFrom purrr transpose
#' @importFrom stats na.omit
#' @export plotUncertainty
#' @author Dongyuan Song, Huy Nyugen

plotUncertainty <- function(ori.tbl,
                            sub.tbl) {

  # Standardize column names
  colnames(ori.tbl) <- c("cell", "ori_pseudotime")
  
  # Combine subsample pseudotimes into a single dataframe
  Merge_pseudotimes <- data.table::rbindlist(sub.tbl, idcol = "subsample")
  
  # Reshape data to wide format
  Merge_pseudotimes <- tidyr::pivot_wider(Merge_pseudotimes, names_from = "subsample", values_from = "pseudotime")
  
  # Rename columns consistently
  colnames(Merge_pseudotimes)[-1] <- paste0("Pseudotime_", seq_len(ncol(Merge_pseudotimes) - 1))
  
  # Rename columns consistently
  colnames(Merge_pseudotimes)[-1] <- paste0("Pseudotime ", seq_len(ncol(Merge_pseudotimes) - 1))
  
  # Merge original pseudotime with subsampled pseudotimes
  Truetimes_pseudotimes <- merge(ori.tbl, Merge_pseudotimes, by = "cell")
  Descending_truetimes_pseudotimes <- dplyr::arrange(Truetimes_pseudotimes, dplyr::desc(ori_pseudotime))
  
  # Extract only pseudotime columns
  pseudotimes_only <- Descending_truetimes_pseudotimes |> dplyr::select(-cell, -ori_pseudotime)
  Descending_pseudotimes_transpose <- purrr::transpose(pseudotimes_only)
  Descending_unlist_pseudotimes <- data.frame(pseudotimes = unlist(Descending_pseudotimes_transpose, use.names = FALSE))
  
  # Replicate cell truetime for each subsample
  cell_truetime_only <- Descending_truetimes_pseudotimes |> dplyr::select(cell, ori_pseudotime)
  multiple_cell_truetime <- cell_truetime_only[rep(seq_len(nrow(cell_truetime_only)), each = length(sub.tbl)), ]
  
  # Combine processed data
  na_plotdata <- cbind(Descending_unlist_pseudotimes, multiple_cell_truetime)
  plotdata <- stats::na.omit(na_plotdata)
  plotdata$cell <- factor(plotdata$cell, levels = unique(plotdata$cell))
  
  # Generate plot
  p <- ggplot(plotdata, aes(x = pseudotimes, y = cell)) +
    stat_density(aes(fill = after_stat(density)), geom = "raster", position = "identity") +
    scale_fill_gradient(low = "white", high = "black") +
    labs(x = "Pseudotimes of subsamples", y = "Cells", fill = "Density") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_blank(),
      aspect.ratio = 1,
      legend.position = "right"
    )
  
  return(p)
}



