#' @title Remove Juveniles
#'
#' @description A function to remove juveniles based on specified threshold.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param char String. Column name of the character to use, typically snout–vent length (SVL).
#' @param threshold Numeric. A value indicating the character value to be considered as adult.
#' @param sum Logical. If TRUE, it prints the number and name of samples removed.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#' 
#' filtered.data <- remove_juv(data, char = "SVL", threshold = 70, sum = TRUE)
#' 
#' filtered.data
#'
#' @return A filtered data frame.
#' @export

remove_juv <- function(data, char = "SVL", threshold = 70, sum = TRUE) {
  ## Access specified column
  col_values <- data[[char]]
  
  ## Identify juveniles based on specified threshold
  juv_rows <- data[col_values < threshold, ]
  
  ## Remove juveniles
  data_clean <- data[col_values >= threshold, ]
  
  ## Print summary
  if (sum) {
    cat("Number of samples removed:", nrow(juv_rows), "\n")
    if (nrow(juv_rows) > 0 && "Pop" %in% names(juv_rows)) {
      cat("Populations of removed samples:\n")
      print(table(juv_rows$Pop))
    }
  }
  
  ## Return the cleaned data
  return(data_clean)
}

