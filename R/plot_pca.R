#' @title Plot Principal Component Analysis
#'
#' @description A function to plot the first two dimensions of Principal Component Analysis (PCA).
#' 
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param fixed.shape A vector specifying the point shape for each population or species. If NULL, random points will be used.
#' @param point.color A vector specifying the point color for each population or species. If NULL, random colors will be used
#' @param ellipse Logical. If TRUE, it adds 95% confidence ellipses. Default TRUE.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' point.shape <- c("Luzon" = 8, "Mindanao" = 11, "Palawan" = 10)
#' point.color <- c("Luzon" = "#000000", "Mindanao" = "#FF7F0E", "Palawan" = "#D62728")
#'
#' pca.plot <- plot_pca(data, point.color = point.color, 
#'                               fixed.shape = point.shape)
#' pca.plot
#'
#' @importFrom stats prcomp
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols rename
#' @importFrom scales breaks_pretty
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 stat_ellipse
#' @return A biplot showing the results of FDA.
#' @export

plot_pca <- function(data, fixed.shape = NULL, point.color = NULL, ellipse = TRUE) {
  ## Ensure 'data' is a data frame and has the necessary columns
  if (!is.data.frame(data)) stop("The input data must be a data frame")
  
  ## Perform PCA
  pca.res <- prcomp(data[2:ncol(data)], center = TRUE, scale. = TRUE)

  ## Get PCA axes values
  pca.pts.load <- as_tibble(pca.res$x) %>% 
    bind_cols(data$Pop) %>%
    rename(Pop = ncol(.))

  ## Prepare PCA data
  pca.data <- subset(pca.pts.load, select = c(ncol(pca.pts.load), 1:2))
  names(pca.data)[1] <- "grp"
  pca.data$grp <- as.factor(pca.data$grp)

  ## Get unique group labels
  pop.lab <- unique(pca.data$grp)
  num.pops <- length(pop.lab)

  ## Set color palette
  if (is.null(point.color)) {
    map.color <- colorRampPalette(brewer.pal(9, "PuBu"))(num.pops)
    names(map.color) <- pop.lab
  } else {
    map.color <- point.color[as.character(pop.lab)]
  }

  ## Set shape palette
  if (is.null(fixed.shape)) {
    set.seed(123)
    point.shape <- sample(1:25, length(pop.lab), replace = TRUE)
    names(point.shape) <- pop.lab
  } else {
    point.shape <- fixed.shape[as.character(pop.lab)]
  }

  ## Plot PCA
  plot_obj <- ggplot(pca.data, aes(x = PC1, y = PC2, color = grp, shape = grp)) +
    geom_point(size = 3) +
    scale_color_manual(values = map.color) + 
    scale_shape_manual(values = point.shape) +
    labs(x = "PCA Component 1", y = "PCA Component 2") +
    theme_classic() +
    theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.background = element_blank(),
      legend.position = "right",
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.text = element_text(size = 15),
      plot.title = element_blank(),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15)
    ) +
    scale_y_continuous(breaks = breaks_pretty()) +
    scale_x_continuous(breaks = breaks_pretty())

  ## Add ellipses
  if (ellipse) {
    plot_obj <- plot_obj +
      stat_ellipse(aes(group = grp), level = 0.95, type = "norm", linetype = 2)
  }

  return(plot_obj)
}

## Declare
utils::globalVariables(c(".", "PC1", "PC2", "grp"))
