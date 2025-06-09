#' @title Principal Component Analysis
#'
#' @description A function to perform Principal Component Analysis (PCA), a  method used to reduce the dimensionality of data by determining
#' the principal components (directions) that maximizes the variance in the dataset. The function uses the 'pincomp' function. It then caculates  
#' the variance explained of the first two components, extracts the variable and point loadings and writes them to csv files.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param prop.var Logical (TRUE or FALSE). If TRUE, it calculates the varaince explained by the first two components.
#' @param var.load Logical (TRUE or FALSE). If TRUE, it extracts and writes the varaible loadings to a csv file.
#' @param pc.scores Logical (TRUE or FALSE). If TRUE, it extracts and writes the PC scores to a csv file.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#' 
#' pca.res <- run_pca(data, prop.var = TRUE, var.load = TRUE, pc.scores = TRUE)
#' 
#' pca.res
#'
#' @importFrom stats prcomp
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols rename
#' @importFrom magrittr %>%
#' @return A list containing the FDA results.
#' @export

run_pca <- function(data, prop.var = TRUE, var.load = TRUE, pc.scores = TRUE) {
  
  ## Perform PCA
  pca.res <- prcomp(data[2:ncol(data)], center = TRUE, scale. = TRUE)
  
  ## Calculate variance explained
  eigs.all <- pca.res$sdev^2
  ax1.var <- eigs.all[1]/sum(eigs.all)*100
  ax1.var <- sprintf(ax1.var, fmt = '%#.2f')
  ax2.var <- eigs.all[2]/sum(eigs.all)*100
  ax2.var <- sprintf(ax2.var, fmt = '%#.2f')

  ## Print proportion of variance explained
  cat(paste0("PC1 variance explained:", ax1.var), "%\n")
  cat(paste0("PC2 variance explained:", ax2.var), "%\n")
  
  ## Extract variable loadings if requested
  pca.var.load <- NULL
  if (var.load) {
    pca.var.load <- as_tibble(pca.res$rotation, rownames = 'char')
    varload_filename <- paste0("pca_varload_", deparse(substitute(data)), ".csv")
    write.csv(pca.var.load, varload_filename)
    # Print only the file write message
    cat(paste0("Variable loadings written to:", varload_filename), "%\n")
  }
  
  ## Extract principal component scores if requested
  pca.pc.scores <- NULL
  if (pc.scores) {
    pca.pc.scores <- as_tibble(pca.res$x) %>% 
      bind_cols(data$Pop) %>%
      rename(Pop = ncol(.))
    pcscores_filename <- paste0("pca_pcscores_", deparse(substitute(data)), ".csv")
    write.csv(pca.pc.scores, pcscores_filename)
    # Print only the file write message
    cat(paste0("Principal component scores written to:", pcscores_filename), "%\n")
  }
  
  ## Return key results
  return(list(
    ax1.var = ax1.var,
    ax2.var = ax2.var,
    pca.var.load = pca.var.load,
    pca.pc.scores = pca.pc.scores
  ))
}


## Declare
utils::globalVariables(".")
