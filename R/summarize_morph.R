#' @title Summarize Morphological Data
#'
#' @description A function to get the summary statistics of all morphological data.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' sum.res <- summarize_morph(data, grp = "Pop")
#'
#' head(sum.res)
#'
#' @importFrom stats median
#' @return A data frame containing summary statistics of all morphological data for every population or species.
#' @export

summarize_morph <- function(data, grp = NULL) {
  ## Check if input is a data frame
  if (!is.data.frame(data)) {
    stop("Error: The input 'data' must be a data frame.")
  }

  if (is.null(grp)) {
    grp <- names(data)[1]
  }

  ## Check if the grp exists in the data frame
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: The column name '", grp, "' does not exist in the data frame.", sep = ""))
  }
  
  ## Check if grp is a factor or character variable
  if (!is.factor(data[[grp]]) && !is.character(data[[grp]])) {
    stop("Error: The grouping variable '", grp, "' must be a factor or character column.")
  }

  groups <- unique(data[[grp]])
  numeric_vars <- names(data)[sapply(data, is.numeric)]
  
  results <- list()
  
  for (group_name in groups) {
    subset_data <- data[data[[grp]] == group_name, ]
    
    for (var in numeric_vars) {
      vals <- subset_data[[var]]
      n_samples <- sum(!is.na(vals))
      
      min_val <- round(min(vals, na.rm = TRUE), 4)
      median_val <- round(median(vals, na.rm = TRUE), 4)
      max_val <- round(max(vals, na.rm = TRUE), 4)
      mean_val <- round(mean(vals, na.rm = TRUE), 4)
      sd_val <- round(sd(vals, na.rm = TRUE), 4)
      Q1 <- round(quantile(vals, 0.25, na.rm = TRUE), 4)
      Q3 <- round(quantile(vals, 0.75, na.rm = TRUE), 4)
      
      stats <- list(
        pop = group_name,
        char = var,
        n = n_samples,
        min = min_val,
        median = median_val,
        max = max_val,
        mean = mean_val,
        sd = sd_val,
        Q1 = Q1,
        Q3 = Q3
      )
      
      results[[length(results) + 1]] <- stats
    }
  }
  
  df <- do.call(rbind, lapply(results, as.data.frame))
  rownames(df) <- NULL

  ## Return data summary
  return(df)
}
