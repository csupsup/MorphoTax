#' @title Summarize Morphological Data
#'
#' @description A function to get the summary statistics of all morphological data.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#' data$Sex <- NULL
#'
#' sum.res <- esult <- summarize_morph(data, grp = "Pop")
#'
#' head(sum.res)
#'
#' @importFrom psych describeBy
#' @return A data frame containing summary statistics of all morphological data for every population or species.
#' @export

summarize_morph <- function(data, grp = "Pop") {
  ## Check if input is a data frame
  if (!is.data.frame(data)) {
    stop("Error: The input 'data' must be a data frame.")
  }
  
  ## Check if the population column exists in the data frame
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: The column name '", grp, "' does not exist in the data frame.", sep = ""))
  }
  
  ## Initialize empty data frame for results
  sum.res <- data.frame()

  ## Loop through columns
  for(col in data[, 2:ncol(data)]) {
    res <- describeBy(col, group = data[[grp]], mat = TRUE)
    sum.res <- rbind(sum.res, res)
  }

  ## Rename character names based on column names
  morph.name <- colnames(data[2:ncol(data)])
  
  ## Get the number of unique populations
  pop <- length(unique(data[[grp]]))
  
  ## Repeat morphological names for each population
  morph.name <- as.data.frame(rep(morph.name, each = pop))
  colnames(morph.name)[1] <- "char"

  ## Add the 'item' column to the summary data frame
  sum.res$item <- morph.name$char

  return(sum.res)
}
