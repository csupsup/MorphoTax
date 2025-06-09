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
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
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
  tukey.summary <- list()

  ## Get names of all morphometric columns (excluding grouping variable)
  morph_names <- setdiff(colnames(data), grp)

  ## Prepare containers for all comparisons and p-values
  all_comparisons <- character()
  tukey_pvals <- list()

  ## Loop through each column
  for (i in seq_along(morph_names)) {
    morph <- morph_names[i]
    col <- data[[morph]]
    
    ## Perform ANOVA
    aov.result <- aov(col ~ data[[grp]], data = data)
    aov.sum <- summary(aov.result)[[1]]
    
    ## Extract ANOVA statistics
    df <- aov.sum$Df[1]
    F.value <- aov.sum$`F value`[1]
    p.value <- aov.sum$`Pr(>F)`[1]
    
    ## Add results to aov.summary
    aov.summary[[i]] <- c(morph, df, F.value, p.value)
    
    ## Perform Tukey HSD test if the p-value is significant
    if (!is.na(p.value) && p.value < 0.05) {
      tukey.result <- TukeyHSD(aov.result)
      tukey.df <- as.data.frame(tukey.result[[1]])
      comparisons <- rownames(tukey.df)
      pvals <- tukey.df$`p adj`

      ## Store comparison names and p-values
      all_comparisons <- unique(c(all_comparisons, comparisons))
      tukey_pvals[[morph]] <- setNames(pvals, comparisons)

      ## Add Comparison and Morph columns for export
      tukey.df$Comparison <- comparisons
      tukey.df$Morph <- morph
      
      ## Add to summary list
      tukey.summary[[length(tukey.summary) + 1]] <- tukey.df
      
      ## Write Tukey result to a CSV file if write.tk is TRUE
      if (write.tk) {
        file_name <- paste0(dir, "tk_", morph, "_results.csv")
        
        ## Ensure the directory exists
        dir.create(dirname(file_name), showWarnings = FALSE, recursive = TRUE)
        
        ## Write the CSV
        write.csv(tukey.df, file = file_name, row.names = FALSE)
      }
    }
  }

  ## Convert ANOVA summary to a data frame
  aov.sum.df <- as.data.frame(do.call(rbind, aov.summary), stringsAsFactors = FALSE)
  colnames(aov.sum.df) <- c("Morph", "DF", "F.value", "p.value")
  aov.sum.df$DF <- as.numeric(aov.sum.df$DF)
  aov.sum.df$F.value <- as.numeric(aov.sum.df$F.value)
  aov.sum.df$p.value <- as.numeric(aov.sum.df$p.value)

  ## Create wide-format matrix: rows = comparisons, columns = morphs, values = p adj
  all_comparisons <- sort(unique(all_comparisons))
  padj_matrix <- matrix(NA, nrow = length(all_comparisons), ncol = length(morph_names))
  colnames(padj_matrix) <- morph_names
  rownames(padj_matrix) <- all_comparisons

  ## Fill in the matrix
  for (morph in names(tukey_pvals)) {
    morph_pvals <- tukey_pvals[[morph]]
    for (comp in names(morph_pvals)) {
      padj_matrix[comp, morph] <- morph_pvals[comp]
    }
  }

  ## Convert to data frame with Comparison as first column
  padj_matrix_df <- as.data.frame(padj_matrix)
  padj_matrix_df <- cbind(Comparison = rownames(padj_matrix_df), padj_matrix_df)
  rownames(padj_matrix_df) <- NULL

  ## Return the list of results
  return(list(
    aov_summary = aov.sum.df,
    tukey_summary = tukey.summary,
    tukey_padj_matrix = padj_matrix_df
  ))
}
