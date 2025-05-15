#' @title Arrange Columns
#'
#' @description A function to reorder columns in morphology data, except for the first one i.e., population or species column.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param col.names A vector containing the names of the columns to reorder. If column names are not included, those will be dropped.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#' cols <- c("SVL", "TD", "OD")
#'
#' data <- arrange_cols(data, col.names = cols)
#' 
#' head(data)
#'
#' @return A data frame containing the reordered columns.
#' @export

arrange_cols <- function(data, col.names) {
  ## Check if all column names in col.names exist in the dataframe
  missing_cols <- setdiff(col.names, colnames(data))
  
  if (length(missing_cols) > 0) {
    ## Print the missing columns
    print(paste("The following columns are missing:", paste(missing_cols, collapse = ", ")))
    stop("Some columns in col.names do not match the columns in the dataframe.")
  }
  
  ## Retain the first column
  first_column <- data[, 1, drop = FALSE]
  
  ## Reorder the remaining columns based on col.names
  remaining_columns <- data[, col.names, drop = FALSE]
  
  ## Combine the first column with the reordered remaining columns
  data <- cbind(first_column, remaining_columns)
  
  ## Return the reordered dataframe
  return(data)
}
