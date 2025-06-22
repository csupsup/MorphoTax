#' @title Plot P-values Bubble Chart
#'
#' @description A function for creating a bubble plot that compares p-values, similar to 'plot_pval', but allows for custom sorting of comparisons.
#' 
#' @param data A data frame with p-values of comparisons.
#' @param dot.val Logical. If TRUE, it adds p-values within the dot.
#' @param sort Logical.If TRUE, sorts the rows and columns.
#' @param custom.sort String specifying the column to use for custom sorting.
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "cleaned.pvals.csv", package = "MorphoTax"))
#'
#' pval.plot <- plot_pval_bubble(data, dot.val = TRUE, sort = TRUE, custom.sort = NULL)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text scale_color_manual scale_size_continuous scale_y_discrete theme_classic theme element_text labs guides guide_legend unit scale_fill_manual
#' @importFrom reshape2 melt
#' @importFrom ggtext element_markdown
#' @importFrom ggnewscale new_scale_fill
#' 
#' @return A bubble plot of p-values
#' @export
#' 

plot_pval_bubble <- function(data, dot.val = TRUE, sort = TRUE, custom.sort = NULL) {

  ## Color-blind–friendly palette (Paul Tol Bright)
  tol_bright <- c(
    "#4477AA", "#66CCEE", "#228833", "#CCBB44",
    "#EE6677", "#AA3377", "#BBBBBB", "#000000",
    "#FFA500", "#00CED1", "#800080", "#DC143C"
  )

  ## Dynamic grouping variable for coloring (based on custom.sort if categorical)
  if (!is.null(custom.sort) && custom.sort %in% colnames(data) &&
      (is.character(data[[custom.sort]]) || is.factor(data[[custom.sort]]))) {
    group_levels <- unique(data[[custom.sort]][order(data[[custom.sort]])])
    n_groups <- length(group_levels)

    if (n_groups <= length(tol_bright)) {
      group_colors <- setNames(tol_bright[1:n_groups], group_levels)
    } else {
      group_colors <- setNames(colorRampPalette(tol_bright)(n_groups), group_levels)
    }
  } else {
    group_levels <- character(0)
    group_colors <- character(0)
  }

  ## Prepare numeric matrix of p-values
  numeric_cols <- sapply(data, is.numeric)
  pval_numeric <- as.matrix(data[, numeric_cols])
  rownames(pval_numeric) <- data$comparisons
  colnames(pval_numeric) <- colnames(data)[numeric_cols]

  ## Reshape to long format
  pval_long <- melt(pval_numeric)
  colnames(pval_long) <- c("Comparison", "Characters", "Pvalue")

  ## Map grouping variable and significance
  if (!is.null(custom.sort) && custom.sort %in% colnames(data)) {
    group_map <- setNames(data[[custom.sort]], data$comparisons)
    pval_long$Group <- group_map[as.character(pval_long$Comparison)]
  } else {
    pval_long$Group <- NA
  }
  pval_long$Significant <- pval_long$Pvalue < 0.05

  ## Sorting logic with flexible custom.sort
  if (!is.null(custom.sort) && custom.sort %in% colnames(data)) {
    sort_col <- data[[custom.sort]]
    if (is.numeric(sort_col)) {
      ordered_comparisons <- data$comparisons[order(sort_col, decreasing = TRUE)]
    } else {
      ordered_comparisons <- data$comparisons[order(as.character(sort_col))]
    }

    sig_counts_cols <- colSums(pval_numeric < 0.05, na.rm = TRUE)
    ordered_characters <- names(sort(sig_counts_cols, decreasing = TRUE))

    pval_long$Comparison <- factor(pval_long$Comparison, levels = ordered_comparisons)
    pval_long$Characters <- factor(pval_long$Characters, levels = ordered_characters)
  } else if (sort) {
    sig_counts_rows <- rowSums(pval_numeric < 0.05, na.rm = TRUE)
    sig_counts_cols <- colSums(pval_numeric < 0.05, na.rm = TRUE)
    ordered_comparisons <- names(sort(sig_counts_rows, decreasing = TRUE))
    ordered_characters <- names(sort(sig_counts_cols, decreasing = TRUE))
    pval_long$Comparison <- factor(pval_long$Comparison, levels = ordered_comparisons)
    pval_long$Characters <- factor(pval_long$Characters, levels = ordered_characters)
  } else {
    pval_long$Comparison <- factor(pval_long$Comparison, levels = data$comparisons)
    pval_long$Characters <- factor(pval_long$Characters, levels = colnames(pval_numeric))
  }

  ## Colored y-axis labels if grouping info available
  if (!is.null(custom.sort) && custom.sort %in% colnames(data) &&
      (is.character(data[[custom.sort]]) || is.factor(data[[custom.sort]]))) {
    comp_group_df <- unique(data[, c("comparisons", custom.sort)])
    y_labels_colored <- sapply(levels(pval_long$Comparison), function(comp) {
      grp <- comp_group_df[[custom.sort]][comp_group_df$comparisons == comp]
      col <- group_colors[grp]
      sprintf("<span style='color:%s'>%s</span>", col, comp)
    })
  } else {
    y_labels_colored <- levels(pval_long$Comparison)
  }

  ## Build base plot
  p <- ggplot(pval_long, aes(x = Characters, y = Comparison)) +
    geom_point(aes(size = 1 - Pvalue, color = Significant)) +
    scale_color_manual(values = c("FALSE" = "#66B2FF", "TRUE" = "#003D73")) +
    scale_size_continuous(range = c(3, 8)) +
    scale_y_discrete(labels = y_labels_colored) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 17),
      axis.text.y = element_markdown(size = 17, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, color = "black"),
      legend.title = element_text(size = 20, face = "bold", color = "black"),
      legend.key.size = unit(1.5, "lines")
    ) +
    labs(
      title = "",
      x = "Characters",
      y = "Comparisons",
      size = "1 - P-value",
      color = "Significant (p < 0.05)",
      fill = custom.sort
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 10)),
      size = guide_legend(order = 2),
      fill = guide_legend(order = 3, override.aes = list(size = 10, alpha = 1))
    )

  ## Optional dot value labels
  if (dot.val) {
    p <- p + geom_text(aes(label = sprintf("%.3f", Pvalue)), size = 3, color = "black")
  }

  ## Add group legend using invisible dummy points if grouping info available
  if (length(group_levels) > 0) {
    p <- p +
      new_scale_fill() +
      geom_point(
        data = data.frame(Group = group_levels),
        aes(x = 0, y = 0, fill = Group),
        shape = 21,
        size = 3,
        color = "black",
        show.legend = TRUE,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(
        name = custom.sort,
        values = group_colors,
        guide = guide_legend(
          override.aes = list(size = 10, alpha = 1)
        )
      )
  }

  return(p)
}

## Declare
utils::globalVariables(c("Group"))
