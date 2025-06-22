#' @title Plot P-values
#'
#' @description A function for creating a dot plot to compare p-values.
#' 
#' @param data A data frame with p-values of comparisons.
#' @param dot.val Logical. If TRUE, it adds p-values within the dot.
#' @param sort Logical.If TRUE, sorts the rows and columns.
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "cleaned.pvals.csv", package = "MorphoTax"))
#'
#' pval.plot <- plot_pval(data, dot.val = TRUE, sort = TRUE)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text scale_color_manual scale_size_continuous
#' @importFrom ggplot2 theme_minimal theme element_text labs
#' @importFrom reshape2 melt
#' @export

plot_pval <- function(data, dot.val = TRUE, sort = TRUE) {
  
  ## Convert to numeric matrix
  pval_numeric <- apply(data, 2, as.numeric)
  rownames(pval_numeric) <- rownames(data)
  colnames(pval_numeric) <- colnames(data)
  
  ## Melt to long format
  pval_long <- melt(pval_numeric, varnames = c("Comparison", "Characters"), value.name = "Pvalue")
  
  ## Flag significant p-values
  pval_long$Significant <- pval_long$Pvalue < 0.05

  ## Sort rows and columns
  if (sort) {
    # Count significant p-values by row and column
    sig_counts_rows <- rowSums(pval_numeric < 0.05, na.rm = TRUE)
    sig_counts_cols <- colSums(pval_numeric < 0.05, na.rm = TRUE)

    # Reorder factors based on significance count
    ordered_comparisons <- names(sort(sig_counts_rows, decreasing = TRUE))
    ordered_characters <- names(sort(sig_counts_cols, decreasing = TRUE))
    
    pval_long$Comparison <- factor(pval_long$Comparison, levels = ordered_comparisons)
    pval_long$Characters <- factor(pval_long$Characters, levels = ordered_characters)
  } else {
    # Keep original order
    pval_long$Comparison <- factor(pval_long$Comparison, levels = rownames(data))
    pval_long$Characters <- factor(pval_long$Characters, levels = colnames(data))
  }

  ## Base plot
  plot_obj <- ggplot(pval_long, aes(x = Characters, y = Comparison)) +
    geom_point(aes(size = 1 - Pvalue, color = Significant)) +
    scale_color_manual(values = c("#66B2FF", "#003D73")) +
    scale_size_continuous(range = c(3, 8)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1),
          axis.text = element_text(size = 15, color = "black"), 
          axis.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    labs(title = "P-values Dot Plot",
         x = "Characters",
         y = "Comparisons",
         size = "1 - P-value",
         color = "Significant (p < 0.05)")

  ## Add labels if dot.val is TRUE
  if (dot.val) {
    plot_obj <- plot_obj + geom_text(aes(label = sprintf("%.3f", Pvalue)), size = 3, color = "black")
  }

  return(plot_obj)
}

utils::globalVariables(c("Characters", "Comparison", "Pvalue", "Significant"))