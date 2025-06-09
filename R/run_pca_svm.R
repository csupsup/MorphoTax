#' @title Support Vector Machine on PCA
#'
#' @description A function to perform support vector machine (SVM) on the first two dimensions of Principal Component Analysis.
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
#' pca.svm <- run_pca_svm(data, point.shape = point.shape, class.color = NULL, pop.order = pop)
#' 
#' pca.svm
#'
#' @importFrom stats prcomp
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom caret trainControl train confusionMatrix
#' @importFrom stats setNames
#' @importFrom ggplot2 geom_tile scale_fill_manual
#' @importFrom scales breaks_pretty
#' @importFrom RColorBrewer brewer.pal
#' @return A plot showing the svm boundary for predicted classes.
#' @export

run_pca_svm <- function(data, point.shape = NULL, class.color = NULL, pop.order = NULL) {

  ## Run PCA and get the first two components
  pca.res <- prcomp(data[2:ncol(data)], center = TRUE, scale. = TRUE)

  # Variance explained
  eigs.all <- pca.res$sdev^2
  ax1.var <- eigs.all[1] / sum(eigs.all) * 100
  ax2.var <- eigs.all[2] / sum(eigs.all) * 100
  cat(sprintf("PC1 variance explained: %.2f%%\n", ax1.var))
  cat(sprintf("PC2 variance explained: %.2f%%\n", ax2.var))

  ## Load PCA scores and add class labels
  pca.pts.load <- as_tibble(pca.res$x[, 1:2]) %>%
    mutate(Pop = data$Pop)

  pca.axes <- pca.pts.load

  ## Order population
  if (!is.null(pop.order)) {
    pca.axes$Pop <- factor(pca.axes$Pop, levels = pop.order)
  }

  ## Set up 10-fold cross-validation
  ctrl <- trainControl(method = "cv", number = 10, classProbs = FALSE)

  ## Train SVM model using caret (linear kernel)
  svm_model <- train(Pop ~ ., data = pca.axes, method = "svmLinear", trControl = ctrl)

  ## Print overall CV results
  print(svm_model)

  ## Use final model to predict and evaluate on full data
  predictions <- predict(svm_model, newdata = pca.axes)
  cm <- confusionMatrix(predictions, pca.axes$Pop)
  print(cm)

  ## Set up plot
  x_range <- seq(min(pca.axes$PC1), max(pca.axes$PC1), length.out = 100)
  y_range <- seq(min(pca.axes$PC2), max(pca.axes$PC2), length.out = 100)
  grid <- expand.grid(PC1 = x_range, PC2 = y_range)

  ## Predict over the PCA grid
  grid$Predicted <- predict(svm_model, newdata = grid)

  ## Default point shapes
  if (is.null(point.shape)) {
    levels_pop <- levels(factor(pca.axes$Pop))
    point.shape <- setNames(1:length(levels_pop), levels_pop)
  }

  ## Default class color
  if (is.null(class.color)) {
    num.pops <- length(unique(grid$Predicted))
    class.color <- colorRampPalette(brewer.pal(9, "PuBu"))(num.pops)
  }

  ## Plot decision boundaries
  plotobj <- ggplot() +
    geom_tile(data = grid, aes(x = PC1, y = PC2, fill = Predicted), alpha = 0.5) +
    geom_point(data = pca.axes, aes(x = PC1, y = PC2, shape = Pop, color = Pop), size = 3, color = "black") +
    scale_shape_manual(values = point.shape) +
    scale_fill_manual(values = class.color) + 
    labs(
      title = "SVM Decision Boundary with PCA Components",
      x = "PCA Component 1",
      y = "PCA Component 2",
      color = "",
      shape = "True Class",
      fill = "Predicted Class"
    ) +
    theme_classic() +
    theme(legend.key = element_blank(),
          legend.title = element_text(size = 15),
          panel.background = element_blank(),
          legend.position = "right",
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.3, "cm"),
          legend.text = element_text(size = 15),
          plot.title = element_blank(),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 15)) +
    scale_y_continuous(breaks = breaks_pretty()) +
    scale_x_continuous(breaks = breaks_pretty()) 
    
  return(plotobj)
}

## Declare
utils::globalVariables(c("PC1", "PC2", "Predicted", "Pop"))
