#' @title Plot Stacked Bar Chart by Characters
#'
#' @description This function counts how many p-values for each character across comparisons are below or above the p-value threshold of 0.05 
#' and then plots a stacked bar chart to visualize the results.
#' 
#' @param data A data frame where the first column is "Comparisons" and the rest are p-values (as character).
#' @param add.num Logical. If TRUE, it adds the total counts on the bars.
#' @param flip Logical. If TRUE, the character labels shifts to the x-axis.
#' @param sort Logical. If TRUE, sorts the character labels based on significance.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "cleaned.pvals.csv", package = "MorphoTax"))
#'
#' pval.stack <- plot_stack_morph(data, add.num = TRUE, flip = FALSE, sort = TRUE)
#' 
#' pval.stack
#' 
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by summarise n %>% mutate
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_fill_manual theme_classic
#' @importFrom ggplot2 theme element_blank element_text unit coord_flip geom_text position_stack
#' @importFrom grid unit
#' 
#' @return A stacked bar plot showing the total count of significant and non-significant characters.
#' @export

plot_stack_morph <- function(data, add.num = TRUE, flip = FALSE, sort = FALSE) {

  ## Convert to long format for character-wise analysis
  data_long <- melt(data, id.vars = "comparisons",
                              variable.name = "character", value.name = "p_value")

  ## Convert p_value to numeric
  data_long$p_value <- as.numeric(data_long$p_value)

  ## Determine significance
  data_long$Significance <- ifelse(data_long$p_value < 0.05, "Significant", "Non-significant")

  ## Count significance by character
  char_counts <- data_long %>%
    group_by(character, Significance) %>%
    summarise(count = n(), .groups = "drop")

  ## Prepare wide format to count significance types per character
  sig_summary <- data_long %>%
    group_by(character) %>%
    summarise(
      less_0.05 = sum(p_value < 0.05, na.rm = TRUE),
      greater_0.05 = sum(p_value >= 0.05, na.rm = TRUE),
      .groups = "drop"
    )

  ## Sort
  if (sort) {
    sig_summary <- sig_summary[order(-sig_summary$less_0.05), ]
  }

  ## Reverse factor levels so that most significant is at the top on y-axis
  char_levels <- rev(sig_summary$character)
  char_counts$character <- factor(char_counts$character, levels = char_levels)

  ## Create plot
  p <- ggplot(char_counts, aes(y = character, x = count, fill = Significance)) +
    geom_bar(stat = "identity") +
    labs(title = "", x = "Number of Comparisons", y = "Character") +
    scale_fill_manual(
      values = c("Significant" = "#003D73", "Non-significant" = "#66B2FF")
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
      plot.title = element_blank(),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 20)
    )

  ## Add labels
  if (add.num) {
    p <- p + geom_text(
      aes(label = count),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 5
    )
  }

  ## Flip axes and tilt labels
  if (flip) {
    p <- p +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }

  return(p)
}

## Declare
utils::globalVariables(c("p_value"))
