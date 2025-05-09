#' @title ANOVA and Tukey HSD
#'
#' @description A function to perform ANOVA on all morphological data across populations or species.
#' Also, if a significant result is detected, it automatically performs Tukey Test.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#' @param write.tk Logical. If TRUE, it writes the results of the Tukey test to csv files.
#' @param dir A directory where the Tukey test results will be stored.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#' data$Sex <- NULL
#'
#' aov.res <- anova_tukey(data, grp = "Pop", write.tk = TRUE, dir = NULL)
#' 
#' aov.res$aov_summary
#' aov.res$tukey_combined
#'
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @importFrom utils write.csv
#'
#' @return A list containing the results of ANOVA and Tukey Test.
#' @export

anova_tukey <- function(data, grp = "Pop", write.tk = FALSE, dir = NULL) {
  ## Ensure the input is a data frame
  if (!is.data.frame(data)) {
    stop("Error: The input 'data' must be a data frame.")
  }

  ## Ensure the population/species column exists in the data
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: The column name '", grp, "' does not exist in the data frame.", sep = ""))
  }
  
  ## Create empty lists to store results
  aov.summary <- list()
  tukey.combined <- list()
  tukey.summary <- list()

  ## Loop through each column
  for (i in 2:ncol(data)) {
    col <- data[, i]
    
    ## Perform ANOVA
    aov.result <- aov(col ~ data[[grp]], data = data)
    aov.sum <- summary(aov.result)[[1]]
    
    ## Extract ANOVA statistics
    df <- aov.sum$Df[1]
    F.value <- aov.sum$`F value`[1]
    p.value <- aov.sum$`Pr(>F)`[1]
    
    ## Add results to aov.summary
    aov.summary[[i - 1]] <- c(colnames(data)[i], df, F.value, p.value)
    
    ## Perform Tukey HSD test if the p-value is significant
    if (p.value < 0.05) {
      tukey.result <- TukeyHSD(aov.result)
      tukey.combined[[i - 1]] <- tukey.result
      
      ## Format Tukey results into a data frame
      tukey.df <- as.data.frame(tukey.result[[1]])
      tukey.df$Comparison <- rownames(tukey.df)
      
      tukey.summary[[i - 1]] <- tukey.df
      
      ## Write Tukey result to a CSV file if write.tk is TRUE
      if (write.tk) {
        col_name <- colnames(data)[i]
        file_name <- paste0(dir, "tk_", col_name, "_results.csv")
        
        ## Ensure the directory exists
        dir.create(dirname(file_name), showWarnings = FALSE, recursive = TRUE)
        
        ## Write the CSV
        write.csv(tukey.df, file = file_name, row.names = FALSE)
      }
    }
  }

  ## Convert ANOVA summary to a data frame
  aov.sum.df <- as.data.frame(do.call(rbind, aov.summary))
  
  ## Set column names for the ANOVA summary
  colnames(aov.sum.df) <- c("Morph", "DF", "F.value", "p.value")
  
  ## Return the list of results
  return(list(
    aov_summary = aov.sum.df,
    tukey_combined = tukey.combined,
    tukey_summary = tukey.summary
  ))
}
