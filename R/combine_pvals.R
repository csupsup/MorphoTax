#' @title Combine P-values
#'
#' @description A function to combine p-values from the functions "aov_tukey" and "kw_dunn".
#' All input CSV files should be located in a folder where each file has a prefix "tk_" for Tukey test results and "dn_" for Dunn test.
#' Each file contains the results of a test on a certain character, for example, "tk_SVL_results.csv."
#'
#' @param dir A directory that contains the results of Tukey and Dunn tests.
#'
#' @examples
#' comb.pvals <- system.file("extdata", "tk_res", package = "MorphoTax")
#'
#' @importFrom dplyr full_join
#' @importFrom utils read.csv
#' @importFrom tools file_path_sans_ext
#' @return A data frame containing the combined p-values.
#' @export

combine_pvals <- function(dir) {

  ## Function to standardize comparisons
  standardize <- function(x) {
    x <- gsub("[-]", " - ", x)
    parts <- strsplit(x, " - ")
    sapply(parts, function(pair) paste(sort(trimws(pair)), collapse = " - "))
  }
  
  ## Get list of CSV files in the directory
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  
  ## Initialize empty list to store data frames
  data_list <- list()
  
  ## Read and standardize each file
  for (file in files) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    
    ## Identify which column has comparisons and which has p-values
    comp_col <- grep("Comparison", names(data), ignore.case = TRUE, value = TRUE)
    pval_col <- grep("p.value|pval|p_value|P.adj|p adj", names(data), ignore.case = TRUE, value = TRUE)
    
    if (length(comp_col) == 0 || length(pval_col) == 0) {
      warning(paste("Skipping file (missing columns):", file))
      next
    }
    
    ## Standardize and rename columns
    data$Normalized_Comparison <- standardize(data[[comp_col]])
    data <- data[, c("Normalized_Comparison", pval_col)]
    
    ## Rename p-value column using file name (without extension)
    pval_name <- file_path_sans_ext(basename(file))
    names(data)[2] <- paste0("pvalue_", pval_name)
    
    data_list[[length(data_list) + 1]] <- data
  }
  
  ## Combine all data frames using full joins
  combined_data <- Reduce(function(x, y) full_join(x, y, by = "Normalized_Comparison"), data_list)
  
  return(combined_data)
}
