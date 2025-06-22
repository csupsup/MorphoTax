#' @title Plot Stacked Bar Chart by Comparison
#'
#' @description This function counts how many p-values for each comparison across are below or above the p-value threshold of 0.05 and
#' then plots a stacked bar chart to visualize the result.
#'
#' @param data A data frame where the first column is "Comparisons" and the rest are p-values (as character).
#' @param add.num Logical. If TRUE, it adds the total counts on the bars.
#' @param flip Logical. If TRUE, the comparison label shifts to the x-axis.
#' @param sort Logical. If TRUE, sorts the comparison based on significane.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "cleaned.pvals.csv", package = "MorphoTax"))
#'
#' pval.stack <- plot_stack_comp(data, add.num = TRUE, flip = FALSE, sort = TRUE)
#' 
#' pval.stack
#' 
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs scale_fill_manual theme_classic theme element_blank
#' @importFrom ggplot2 element_text unit coord_flip
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom stats na.omit
#' @importFrom grid unit
#' 
#' @return A stacked bar plot showing the total count of significant and non-significant characters.
#' @export

plot_stack_comp <- function(data, add.num = TRUE, flip = FALSE, sort = FALSE) {

  ## Convert p-values to numeric
  data_numeric <- data
  data_numeric[, -1] <- lapply(data[, -1], as.numeric)

  ## Count significant and non-significant p-values for each comparison
  pval_counts <- data.frame(
    comparisons = data_numeric$comparisons,
    less_0.05 = apply(data_numeric[, -1], 1, function(x) sum(x < 0.05, na.rm = TRUE)),
    greater_0.05 = apply(data_numeric[, -1], 1, function(x) sum(x >= 0.05, na.rm = TRUE))
  )

  ## Sort by number of significant p-values
  if (sort) {
    pval_counts <- pval_counts[order(-pval_counts$less_0.05), ]
  }

  ## Reverse factor levels so comparisons with most significance are on top of y-axis
  comp_levels <- rev(pval_counts$comparisons)
  pval_counts$comparisons <- factor(pval_counts$comparisons, levels = comp_levels)

  ## Reshape to long format for ggplot
  pval_long <- melt(
    pval_counts[, c("comparisons", "less_0.05", "greater_0.05")],
    id.vars = "comparisons",
    variable.name = "Significance", value.name = "count"
  )

  ## Non-significant at bottom, Significant on top
  pval_long$Significance <- factor(pval_long$Significance, levels = c("greater_0.05", "less_0.05"))

  ## Create plot (stacking happens along x-axis even if y-axis is categorical)
  p <- ggplot(pval_long, aes(y = comparisons, x = count, fill = Significance)) +
    geom_bar(stat = "identity") +
    labs(x = "Number of Characters", y = "Comparison") +
    scale_fill_manual(
      values = c("greater_0.05" = "#66B2FF", "less_0.05" = "#003D73"),
      labels = c("Non-significant", "Significant")
    ) +
    theme_classic() +
    theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.background = element_blank(),
      legend.position = "top",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 20)
    )

  ## Add text labels to bars
  if (add.num) {
    p <- p + geom_text(
      aes(label = count),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 5
    )
  }

  ## Flip axes
  if (flip) {
    p <- p +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }

  return(p)
}

## Declare
utils::globalVariables(c("comparisons", "count", "Significance"))