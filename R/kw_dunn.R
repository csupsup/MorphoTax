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
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
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
  dunn.combined <- list()
  dunn.summary <- list()

  ## Loop through each column, except the first one.
  for (i in 2:ncol(data)) {
    col <- data[, i]
    
    ## Perform Kruskal-Wallis test
    kruskal.result <- kruskal.test(col ~ data[[grp]], data = data)
    kruskal.p.value <- kruskal.result$p.value
    
    ## Store Kruskal-Wallis results
    kruskal.summary[[i - 1]] <- c(colnames(data)[i], kruskal.result$statistic, kruskal.p.value)
    
    ## Perform Dunn's Test if Kruskal-Wallis p-value is significant
    if (kruskal.p.value < 0.05) {
      dunn.result <- dunnTest(col ~ data[[grp]], data = data, method = "bonferroni")
      dunn.combined[[i - 1]] <- dunn.result
      
      dunn.df <- as.data.frame(dunn.result$res)
      dunn.summary[[i - 1]] <- dunn.df
      
      ## Write Dunn's result to a csv file if write.dunn is TRUE
      if (write.dunn) {
        col.name <- colnames(data)[i]
        file.name <- paste0(dir, "dn_", col.name, "_results.csv")
        
        ## Ensure the directory exists
        dir.create(dirname(file.name), showWarnings = FALSE, recursive = TRUE)
        
        ## Write the csv
        write.csv(dunn.df, file = file.name, row.names = FALSE)
      }
    }
  }

  ## Convert Kruskal-Wallis summary to a data frame
  kruskal.sum <- as.data.frame(do.call(rbind, kruskal.summary))
  
  ## Add column names to the Kruskal-Wallis result
  colnames(kruskal.sum) <- c("Morph", "KW Statistic", "p.value")
  
  ## Return the list of results
  return(list(
    kruskal_summary = kruskal.sum,
    dunn_combined = dunn.combined,
    dunn_summary = dunn.summary
  ))
}
