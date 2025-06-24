#' @title Support Vector Machine on FDA
#'
#' @description A function to perform support vector machine (SVM) on the first two dimensions of Flexible Discriminant Analysis.
#' SVM is a classification algorithm that partitions samples into classes by searching the optimal hyperplane margin 
#' that maximizes the boundary among classes (in this case, putative population or species) in two or more dimensions.
#' The function uses the 'svm' function from the package 'e1071', implementing the 'svmLinear'method and a 10-fold cross-validation. 
#' The function ouputs a plot showing the svm boundary for predicted classes.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param point.shape A vector specifying the point shapes for populations. If NULL, it assigns random point shapes.
#' @param class.color A vector specifying the point and class colors. If NULL, it assigns random point shapes and colors.
#' @param pop.order A vector specifying the order of populations to be shown in the plot legend.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' point.shape <- c("Luzon" = 8, "Mindanao" = 11, "Palawan" = 10)
#' point.color <- c("Luzon" = "#000000", "Mindanao" = "#FF7F0E", "Palawan" = "#D62728")
#' pop <- c("Palawan", "Mindanao", "Luzon")
#' 
#' fda.svm <- run_fda_svm(data, point.shape = point.shape, class.color = NULL, pop.order = pop)
#' fda.svm
#'
#' @importFrom caret trainControl train confusionMatrix
#' @importFrom stats setNames
#' @importFrom ggplot2 geom_tile scale_fill_manual
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales breaks_pretty
#' @return A plot showing the svm boundary for predicted classes.
#' @export

run_fda_svm <- function(data, point.shape = NULL, class.color = NULL, pop.order = pop) {

  ## Run FDA and extract first two dimensions
  fda.res <- fda(Pop ~ ., data = data)
  fda_axes <- as.data.frame(fda.res$fit$fitted.values[, 1:2])
  fda_axes <- cbind(Pop = data$Pop, fda_axes)

  ## Factor ordering if provided
  if (!is.null(pop.order)) {
    fda_axes$Pop <- factor(fda_axes$Pop, levels = pop.order)
  }

  ## Rename dimensions
  colnames(fda_axes)[2:3] <- c("FDA1", "FDA2")

  ## Set up 10-fold cross-validation
  ctrl <- trainControl(method = "cv", number = 10, classProbs = FALSE)

  ## Train SVM on FDA-reduced data
  svm_model <- train(Pop ~ ., data = fda_axes, method = "svmLinear", trControl = ctrl)

  ## Print cross-validated performance
  print(svm_model)

  ## Evaluate predictions on full dataset for confusion matrix
  predictions <- predict(svm_model, newdata = fda_axes)
  cm <- confusionMatrix(predictions, fda_axes$Pop)
  print(cm)

  ## Create grid for decision boundary
  x_range <- seq(min(fda_axes$FDA1), max(fda_axes$FDA1), length.out = 100)
  y_range <- seq(min(fda_axes$FDA2), max(fda_axes$FDA2), length.out = 100)
  grid <- expand.grid(FDA1 = x_range, FDA2 = y_range)
  grid$Predicted <- predict(svm_model, newdata = grid)

  ## Set default point shapes
  if (is.null(point.shape)) {
    levels_pop <- levels(factor(fda_axes$Pop))
    point.shape <- setNames(1:length(levels_pop), levels_pop)
  }

  ## Set default class color
  if (is.null(class.color)) {
    num.pops <- length(unique(grid$Predicted))
    class.color <- colorRampPalette(brewer.pal(9, "PuBu"))(num.pops)
  }

  ## Plot decision boundaries on FDA dimensions
  plotobj <- ggplot() +
    geom_tile(data = grid, aes(x = FDA1, y = FDA2, fill = Predicted), alpha = 0.5) +
    geom_point(data = fda_axes, aes(x = FDA1, y = FDA2, shape = Pop), size = 3, color = "black") +
    scale_shape_manual(values = point.shape) +
    scale_fill_manual(values = class.color) +
    labs(
      title = "SVM Decision Boundary with FDA Components",
      x = "FDA Dimension 1",
      y = "FDA Dimension 2"
    ) +
    guides(
      shape = guide_legend(order = 1, title = "True Class"),
      fill = guide_legend(order = 2, title = "Predicted Class")
    ) +
    theme_classic() +
    theme(
      legend.key = element_blank(),
      legend.title = element_text(size = 15),
      panel.background = element_blank(),
      legend.position = "right",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.text = element_text(size = 15),
      plot.title = element_blank(),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15)
    ) +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    scale_x_continuous(breaks = scales::breaks_pretty()) 

  return(plotobj)
}

## Declare
utils::globalVariables(c("pop", "FDA1", "FDA2", "Predicted", "Pop"))
