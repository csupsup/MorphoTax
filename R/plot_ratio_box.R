#' @title Ratio Boxplot
#'
#' @description A function to create boxplot of character ratio.
#' 
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param char.as.num String. Column name for character representing the numerator.
#' @param char.as.den String. Column name for character representing the denominator.
#' @param grp  String. Column name for population or species.
#' @param y.title Custom label for y-axis.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' 
#' box <- plot_ratio_box(data, char.as.num = "OD", char.as.den = "SVL", grp = "Pop", 
#'                y.title = "Orbit Diameter/\nSnout–Vent Length")
#' box
#'
#' @importFrom ggplot2 geom_boxplot margin
#' @return A boxplot of ratio.
#' @export

plot_ratio_box <- function(data, char.as.num, char.as.den, grp, y.title = NULL) {
  
  ## Check if char.as.num, char.as.den, and grp are present in the data frame
  missing_cols <- setdiff(c(char.as.num, char.as.den, grp), colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following variables are missing from the data frame:", paste(missing_cols, collapse = ", ")))
  }
  
  ## Calculate the ratio of char.as.num and char.as.den
  data$ratio <- data[[char.as.num]] / data[[char.as.den]]
  
  ## Create the boxplot
  boxplot <- ggplot(data, aes(x = .data[[grp]], y = .data$ratio)) +
    theme_classic() +
    geom_boxplot(fill = "#5B8CB6") +
    labs(title = "", 
         x = "", 
         y = y.title %||% paste(char.as.num, "/", "\n", char.as.den)) + 
    theme(axis.text = element_text(size = 15, color = "black"), 
          axis.title = element_text(size = 15, face = "bold"), 
          plot.title = element_blank(),
          legend.background = element_blank(),
          axis.text.x = element_text(angle = 50, hjust = 1), 
          axis.text.y = element_text(size = 12),
          plot.margin = margin(10, 0, 0, 0),
          axis.title.x = element_blank())
  
  ## Return the plot
  return(boxplot)
}

## Declare
utils::globalVariables(c(".data"))