#' @title Remove Non-Numeric Samples
#'
#' @description A function to remove samples (rows) with non-numeric or empty cells.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param cols.range Vector. Column range to check for non-numeric or empty cells. Input can be "3:6" or "3:ncol(data)". Default 3:ncol(data).
#' @param sum Logical. If TRUE, it prints the number and name of samples removed.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' 
#' filtered.data <- remove_non_num(data, cols.range = 3:ncol(data), sum = TRUE)
#' 
#' filtered.data
#'
#' @return A filtered data frame.
#' @export

remove_non_num <- function(data, cols.range = 3:ncol(data), sum = TRUE) {
  ## Ensure column inputs is correct
  if (is.character(cols.range)) {
    cols.range <- as.numeric(cols.range)
  }
  
  ## Identify non-numeric or samples with empty cells
  non_numeric_or_empty_rows <- apply(data[, cols.range], 1, function(row) {
    any(sapply(row, function(x) {
      is.character(x) && (is.na(suppressWarnings(as.numeric(x))) || x == "")
    }))
  })
  
  ## Extract sample names
  removed_samples <- data[non_numeric_or_empty_rows, 1]  
  
  ## Record the number of samples removed
  num_removed <- sum(non_numeric_or_empty_rows)
  
  ## Remove those samples
  clean_data <- data[!non_numeric_or_empty_rows, ]
  
  ## Convert specified columns to numeric
  clean_data[, cols.range] <- lapply(clean_data[, cols.range], function(col) {
    as.numeric(col)
  })
  
  ## Print summary
  if (sum) {
    cat("Number of samples removed:", num_removed, "\n")
    cat("Remaining number of samples:", nrow(clean_data), "\n")
    cat("Removed sample names:\n")
    print(removed_samples)
  }
  
  ## Return the cleaned data
  return(clean_data)
}
