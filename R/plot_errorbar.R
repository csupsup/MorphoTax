#' @title Plot Error Bar
#'
#' @description A function for creating an error bar plot and test differences in characters between 
#' sexes within each population using the T-test and Wilcoxon test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param char String. Column name for morphological character to test.
#' @param grp String. Column name for population or species.
#' @param test String. Specify which test to use ("t-test", "wilcox", "none").
#' @param pop.order A vector specifying the order of populations to be shown in the plot legend.
#' @param asterisk Logical. If TRUE, it adds an asterisk to the plot to indicate a significant result.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' pop <- c("Luzon", "Palawan", "Mindanao")
#' 
#' svl.bar <- plot_errorbar(data, char = "SVL", sex = "Sex", 
#'              grp = "Pop", test = "t-test", pop.order = pop, asterisk = TRUE)
#' 
#' svl.bar
#'
#' @importFrom ggplot2 ggplot aes aes_string geom_errorbar geom_line geom_point 
#' @importFrom ggplot2 scale_color_manual scale_shape_manual ylab xlab theme theme_classic
#' @importFrom ggplot2 element_text element_blank geom_text
#' @importFrom Rmisc summarySE
#' @importFrom stats t.test wilcox.test
#' @importFrom rlang sym
#' 
#' @export

plot_errorbar <- function(data, char = "SVL", sex = "Sex", grp = "Pop", 
                          test = c("t-test", "wilcox", "none"), 
                          pop.order = NULL, 
                          asterisk = TRUE) {
  test <- match.arg(test)
  
  ## Check required columns exist
  required_cols <- c(char, sex, grp)
  if (!all(required_cols %in% names(data))) {
    stop("One or more specified columns (char, sex, grp) do not exist in the data.")
  }
  
  ## Summary statistics with standard error
  df.sum <- summarySE(data, measurevar = char, groupvars = c(sex, grp))
  df.sum[[grp]] <- factor(df.sum[[grp]])
  df.sum[[sex]] <- factor(df.sum[[sex]])
  
  ## Calculate errorbar limits
  df.sum$ymin <- df.sum[[char]] - df.sum$se
  df.sum$ymax <- df.sum[[char]] + df.sum$se
  
  sig_annotations <- data.frame()

  if (asterisk && test != "none") {
    group_levels <- unique(data[[grp]])
    
    for (g in group_levels) {
      sub_data <- data[data[[grp]] == g, ]
      if (length(unique(sub_data[[sex]])) == 2) {
        
        ## Perform the chosen test
        p_val <- switch(
          test,
          "t-test" = t.test(sub_data[[char]] ~ sub_data[[sex]])$p.value,
          "wilcox" = wilcox.test(sub_data[[char]] ~ sub_data[[sex]])$p.value
        )
        
        if (p_val < 0.05) {
          pos <- df.sum[df.sum[[grp]] == g, ]
          max_y <- max(pos[[char]] + pos$se)
          offset <- 0.05 * max_y
          
          sig_annotations <- rbind(sig_annotations, data.frame(
            grp_val = g,
            y = max_y + offset,
            label = "*"
          ))
        }
      }
    }
  }
  
  ## Prepare significance annotation dataframe for plotting
  if (nrow(sig_annotations) > 0) {
    sig_annotations[[grp]] <- sig_annotations$grp_val
  }
  
  ## Order population
  if (!is.null(pop.order)) {
    df.sum[[grp]] <- factor(df.sum[[grp]], levels = pop.order)
    if (nrow(sig_annotations) > 0) {
      sig_annotations[[grp]] <- factor(sig_annotations[[grp]], levels = pop.order)
    }
  }

  ## Plot
  plot_obj <- ggplot(df.sum, aes(x = !!sym(grp), y = !!sym(char), colour = !!sym(sex), group = !!sym(sex))) +
    theme_classic() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1) +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = !!sym(sex)), size = 3) +
    scale_color_manual(name = "", labels = c("Female", "Male"), values = c("#66B2FF", "#003D73")) +
    scale_shape_manual(name = "", labels = c("Female", "Male"), values = c(15, 16)) +
    ylab(char) +
    xlab("") +
    theme(
      axis.text = element_text(size = 15, color = "black"), 
      axis.title = element_text(size = 15, face = "bold"), 
      legend.text = element_text(size = 15),
      legend.position = c(0.80, 0.95),
      legend.background = element_blank(),
      axis.text.x = element_text(angle = 50, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_blank()
    )

  ## Add significance asterisks if enabled and present
  if (asterisk && nrow(sig_annotations) > 0) {
    plot_obj <- plot_obj + geom_text(data = sig_annotations, 
                                     aes(x = !!sym(grp), y = y, label = label), 
                                     inherit.aes = FALSE, size = 6, vjust = 0)
  }

  return(plot_obj)
}

## Declare
utils::globalVariables(c("ymin", "ymax" ,"y", "label"))