#' @title Plot Flexible Discriminant Analysis
#'
#' @description A function to plot the first two dimensions of Flexible Discriminant Analaysis (FDA).
#' 
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param fixed.shape A vector specifying the point shape for each population or species. If NULL, random points will be used.
#' @param point.color A vector specifying the point color for each population or species. If NULL, random colors will be used
#' @param split Numeric. Proportion of data to be assigned as training data. If NULL, all data will be used as training set.
#' @param accu Logical. If TRUE, it prints the accuracy of classification based on confusion matrix.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' point.shape <- c("Luzon" = 8, "Mindanao" = 11, "Palawan" = 10)
#' point.color <- c("Luzon" = "#000000", "Mindanao" = "#FF7F0E", "Palawan" = "#D62728")
#'
#' fda.plot <- plot_fda(data, point.color = point.color, 
#'                                fixed.shape = point.shape, split = 0.7, accu = TRUE)
#' fda.plot
#'
#' @importFrom mda fda
#' @importFrom stats predict
#' @importFrom car leveneTest
#' @importFrom stats as.formula
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot geom_point aes scale_color_manual
#'   scale_shape_manual labs theme_classic theme element_blank
#'   element_text scale_y_continuous scale_x_continuous
#' @importFrom grid unit
#' @return A biplot showing the results of FDA.
#' @export

plot_fda <- function(data, fixed.shape = NULL, point.color = NULL, split = NULL, accu = TRUE) {
  if (!is.data.frame(data)) stop("The input data must be a data frame")
  if (!"Pop" %in% names(data)) stop("The data must contain a 'Pop' column")

  ## Data split  
  if (is.null(split)) {
    train_data <- data
    test_data <- data
  } else {
    shuffled_index <- sample(1:nrow(data))
    train_size <- floor(split * nrow(data))
    train_index <- shuffled_index[1:train_size]
    test_index <- shuffled_index[(train_size + 1):nrow(data)]
    train_data <- data[train_index, ]
    test_data <- data[test_index, ]
  }

  ## Perform FDA on training data 
  fda.res <- fda(Pop ~ ., data = train_data)

  ## Get FDA axes (project both train and test data)
  fda.axes <- predict(fda.res, newdata = test_data, type = "var")
  
  ## Convert fda.axes to a data frame
  fda.axes2 <- as.data.frame(fda.axes)
  
  ## Check the structure of fda.axes2
  if (ncol(fda.axes2) >= 2) {
    colnames(fda.axes2)[1:2] <- c("V1", "V2")
  } else {
    stop("FDA axes must contain at least two dimensions")
  }

  grp <- test_data$Pop
  fda.df <- cbind(grp, fda.axes2)
  fda.df <- subset(fda.df, select = c(1:3))
  fda.df$grp <- as.factor(fda.df$grp)

  ## Color and shape setup 
  pop.lab <- unique(fda.df$grp)
  num.pops <- length(pop.lab)

  if (is.null(point.color)) {
    map.color <- colorRampPalette(brewer.pal(9, "PuBu"))(num.pops)
  } else {
    map.color <- point.color[as.character(pop.lab)]
  }

  if (is.null(fixed.shape)) {
    set.seed(123)
    point.shape <- sample(1:25, length(pop.lab), replace = TRUE)
  } else {
    point.shape <- fixed.shape[as.character(pop.lab)]
  }

  ## Plot
  plot_obj <- ggplot() +
    geom_point(data = fda.df, aes(x = V1, y = V2, color = grp, shape = grp), size = 3) +
    scale_color_manual(values = map.color) + 
    scale_shape_manual(values = point.shape) +
    labs(x = "FDA Dimension 1", y = "FDA Dimension 2") +
    theme_classic() +
    theme(legend.key = element_blank(), 
          legend.title = element_blank(), 
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

  ## Accuracy
  if (accu) {
    if (is.null(split)) {
      pred <- predict(fda.res, newdata = train_data)
      actual <- train_data$Pop
      label <- "Training"
    } else {
      pred <- predict(fda.res, newdata = test_data)
      actual <- test_data$Pop
      label <- "Test"
    }

    cm <- table(actual, pred)
    accuracy <- sum(diag(cm)) / sum(cm) * 100
    cat(paste0("FDA Classification Accuracy on ", label, " Data: ", round(accuracy, 2), "%\n"))
  }

  return(plot_obj)
}

## Declare
utils::globalVariables(c("V1", "V2"))
