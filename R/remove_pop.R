#' @title Remove Populations
#'
#' @description A function to remove populations with small samples or individuals.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#' @param threshold Numeric. Number of samples required for each population to be retained.
#' Populations with fewer than this number will be removed. Default is 4.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#' 
#' filtered.data <- remove_pop(data, grp = "Pop", threshold = 4)
#' 
#' filtered.data
#'
#' @return A filtered data frame.
#' @export

remove_pop <- function(data, grp = "Pop", threshold = 4) {
  ## Check if 'grp' exists in the data frame
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "is not in the data frame."))
  }
  
  ## Check if 'threshold' is a valid numeric value and >= 4
  if (is.null(threshold) || !is.numeric(threshold) || threshold < 4) {
    stop("Error: 'threshold' must be a numeric value and >= 4.")
  }
  
  ## Identify populations with fewer than 'threshold' individuals
  pop_counts <- table(data[[grp]])
  populations_to_remove <- names(pop_counts[pop_counts < threshold])
  
  ## Print which populations are being removed
  if (length(populations_to_remove) > 0) {
    cat("Removing populations with fewer than", threshold, "individuals:", 
        paste(populations_to_remove, collapse = ", "), "\n")
  } else {
    cat("No populations were removed. All populations have at least", threshold, "individuals.\n")
  }
  
  ## Remove populations with < 'threshold' individuals
  filtered.data <- data[data[[grp]] %in% 
                        names(pop_counts[pop_counts >= threshold]),]
  
  return(filtered.data)
}
