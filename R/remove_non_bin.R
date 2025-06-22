#' @title Remove Samples with Non-binary Values
#'
#' @description A function to remove samples with non-binary values.
#' 
#' @param data A data frame to filter.
#' @param grp String. Name of group or population to report number of removed samples.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' 
#' filtered.data <- remove_non_bin(data, grp = "Pop")
#' 
#' head(filtered.data)
#' 
#' @return A filtered data frame with binary data.
#' @export

remove_non_bin <- function(data, grp = NULL) {
  ## Check if input is a data frame
  if (!is.data.frame(data)) {
    stop("Error: Input must be a data frame.")
  }
  
  ## Identify numeric columns
  num_cols <- sapply(data, is.numeric)
  
  ## Logical vector: TRUE if row has any numeric value >1 or <0
  rows_to_remove <- apply(data[, num_cols, drop = FALSE], 1, function(x) any(x > 1 | x < 0))
  
  if (!is.null(grp)) {
    ## Check group column 
    if (!grp %in% names(data)) {
      stop(paste0("Error: Column '", grp, "' not found in data frame."))
    }
    
    ## Count removals per group
    removed_counts <- table(data[[grp]][rows_to_remove])
    ## Print counts
    all_groups <- unique(data[[grp]])
    for (g in all_groups) {
      count <- ifelse(g %in% names(removed_counts), removed_counts[[g]], 0)
      message(sprintf("Removed %d rows from pop '%s'", count, g))
    }
  }
  
  ## Filter data
  filtered_data <- data[!rows_to_remove, ]
  
  ## Return filtered data
  return(filtered_data)
}
