#' @title Kruskal-Wallis and Dunn Test
#'
#' @description A function to perform Kruskal-Wallis test on all morphological data.
#' Also, if a significant result is detected, it automatically performs Dunn Test with a Bonferroni adjustment for p-values.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#' @param write.dunn Logical(TRUE or FALSE). If TRUE, it writes the results of the Dunn test to csv files.
#' @param dir A directory where the Dunn test results will be stored.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' kw.res <- kw_dunn(data, grp = "Pop", write.dunn = TRUE, dir = NULL)
#' 
#' ## access results
#' kw.res$kruskal_summary
#' kw.res$dunn_combined
#'
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom stats kruskal.test
#' @importFrom FSA dunnTest
#' @return A list containing the results of Kruskal-Wallis and Dunn Test.
#' @export

kw_dunn <- function(data, grp = "Pop", write.dunn = TRUE, dir = "dn_female/dn_") {
  ## Ensure the input is a data frame
  if (!is.data.frame(data)) {
    stop("Error: The input 'data' must be a data frame.")
  }
  
  ## Ensure the population/species column exists in the data
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: The column name '", grp, "' does not exist in the data frame.", sep = ""))
  }
  
  data[[grp]] <- as.factor(data[[grp]])

  ## Create empty lists to store results
  kruskal.summary <- list()
  dunn.summary <- list()

  ## Prepare containers for comparison names and p-values
  morph_names <- setdiff(colnames(data), grp)
  all_comparisons <- character()
  dunn_pvals <- list()

  ## Loop through each column
  for (i in seq_along(morph_names)) {
    morph <- morph_names[i]
    col <- data[[morph]]
    
    ## Perform Kruskal-Wallis test
    kruskal.result <- kruskal.test(col ~ data[[grp]], data = data)
    kruskal.p.value <- kruskal.result$p.value
    
    ## Store Kruskal-Wallis results
    kruskal.summary[[i]] <- c(morph, kruskal.result$statistic, kruskal.p.value)
    
    ## Perform Dunn's Test if Kruskal-Wallis p-value is significant
    if (!is.na(kruskal.p.value) && kruskal.p.value < 0.05) {
      dunn.result <- dunnTest(col ~ data[[grp]], data = data, method = "bonferroni")
      dunn.df <- as.data.frame(dunn.result$res)
      
      ## Extract comparison names and p-values
      comparisons <- dunn.df$Comparison
      pvals <- dunn.df$P.adj
      
      ## Store in lookup list for padj matrix
      all_comparisons <- unique(c(all_comparisons, comparisons))
      dunn_pvals[[morph]] <- setNames(pvals, comparisons)
      
      ## Add Morph column for record keeping
      dunn.df$Morph <- morph
      dunn.summary[[length(dunn.summary) + 1]] <- dunn.df

      ## Write Dunn's result to a csv file if write.dunn is TRUE
      if (write.dunn) {
        file.name <- paste0(dir, "dn_", morph, "_results.csv")
        
        ## Ensure the directory exists
        dir.create(dirname(file.name), showWarnings = FALSE, recursive = TRUE)
        
        ## Write the csv
        write.csv(dunn.df, file = file.name, row.names = FALSE)
      }
    }
  }

  ## Convert Kruskal-Wallis summary to a data frame
  kruskal.sum <- as.data.frame(do.call(rbind, kruskal.summary), stringsAsFactors = FALSE)
  colnames(kruskal.sum) <- c("char", "KW Statistic", "p.value")
  kruskal.sum$`KW Statistic` <- as.numeric(kruskal.sum$`KW Statistic`)
  kruskal.sum$p.value <- as.numeric(kruskal.sum$p.value)

  ## Create p.adj matrix: rows = comparisons, columns = morphs
  all_comparisons <- sort(unique(all_comparisons))
  padj_matrix <- matrix(NA, nrow = length(all_comparisons), ncol = length(morph_names))
  colnames(padj_matrix) <- morph_names
  rownames(padj_matrix) <- all_comparisons

  for (morph in names(dunn_pvals)) {
    morph_pvals <- dunn_pvals[[morph]]
    for (comp in names(morph_pvals)) {
      padj_matrix[comp, morph] <- morph_pvals[comp]
    }
  }

  ## Convert to data frame with Comparison column
  padj_matrix_df <- as.data.frame(padj_matrix)
  padj_matrix_df <- cbind(Comparison = rownames(padj_matrix_df), padj_matrix_df)
  rownames(padj_matrix_df) <- NULL

  ## Return the list of results
  return(list(
    kruskal_summary = kruskal.sum,
    dunn_summary = dunn.summary,
    dunn_padj_matrix = padj_matrix_df
  ))
}
