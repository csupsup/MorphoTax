#' @title Remove Outliers
#'
#' @description A function to remove outliers based on interquartile range method.
#' It uses only one character as a reference in removing samples.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param char String. Column name for character to use for removing outliers.
#' @param grp String. Column name for population or species. If NULL, the entire dataset will be considered as one population.
#' @param q1 Numeric. Value for lower quartile. Default is set to 0.25.
#' @param q3 Numeric. Value for upper quartile. Default is set to 0.75.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#' 
#' filtered.data <- remove_outliers(data, char = "SVL", grp = "Pop", q1 = 0.25, q3 = 0.75)
#' 
#' filtered.data
#'
#' @importFrom stats quantile
#' @return A filtered data frame.
#' @export

remove_outliers <- function(data, char = "SVL", grp = NULL, q1 = 0.25, q3 = 0.75) {
  if (is.null(grp)) {
    ## If no group is provided, treat all data as one group
    Q1_value <- quantile(data[[char]], q1, na.rm = TRUE)
    Q3_value <- quantile(data[[char]], q3, na.rm = TRUE)
    IQR_value <- Q3_value - Q1_value

    ## Calculate the lower and upper bounds for the entire dataset
    lower_bound <- Q1_value - 1.5 * IQR_value
    upper_bound <- Q3_value + 1.5 * IQR_value

    ## Identify which rows are outliers
    is_outlier <- data[[char]] < lower_bound | data[[char]] > upper_bound

    ## Calculate number of samples removed
    num_removed <- sum(is_outlier, na.rm = TRUE)

    ## Message indicating outlier removal across the entire dataset
    message("Outlier removal applied across the entire dataset. ", num_removed, " samples removed.")
  } else {
    ## If a grouping chariable is provided, calculate bounds for each group
    Q1_values <- tapply(data[[char]], data[[grp]], function(x) quantile(x, q1, na.rm = TRUE))
    Q3_values <- tapply(data[[char]], data[[grp]], function(x) quantile(x, q3, na.rm = TRUE))
    IQR_values <- Q3_values - Q1_values

    ## Calculate the lower and upper bounds for each group
    lower_bound <- Q1_values - 1.5 * IQR_values
    upper_bound <- Q3_values + 1.5 * IQR_values

    ## Identify which rows are outliers
    is_outlier <- data[[char]] < lower_bound[data[[grp]]] | data[[char]] > upper_bound[data[[grp]]]

    ## Calculate the number of samples removed per group
    num_removed_per_group <- tapply(is_outlier, data[[grp]], sum, na.rm = TRUE)

    ## Print number of samples removed per group
    message("Outlier removal applied per ", grp,".")
    for (group in names(num_removed_per_group)) {
      message("Group '", group, "' removed ", num_removed_per_group[group], " samples.")
    }
  }

  ## Filter out the outliers and return the cleaned dataset
  data_clean <- data[!is_outlier, ]

  return(data_clean)
}
