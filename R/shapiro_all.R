#' @title Shapiro-Wilk Test
#'
#' @description A function to test normality of all morphological data using Shapiro-Wilk Test.
#' 
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#' 
#' shapiro.res <- shapiro_all(data, grp = "Pop")
#'
#' shapiro.res
#'
#' @importFrom stats shapiro.test
#' @return A data frame summarizing the results of Shapiro-Wilk Test.
#' @export

shapiro_all <- function(data, grp = "Pop") {
  ## Check if the specified column exists in the data frame
  if (!is.null(grp) && !(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the data frame."))
  }

  ## Initialize a list to store the results
  results <- list()

  ## Check if grouping is specified (i.e., grp is not NULL)
  if (is.null(grp)) {
    ## Perform Shapiro-Wilk test on each numeric column across the entire data
    for (char in names(data)) {
      if (is.numeric(data[[char]])) {
        if (length(unique(data[[char]])) > 1) {
          shapiro.res <- shapiro.test(data[[char]])  
          p.value <- shapiro.res$p.value  
          significance <- ifelse(p.value < 0.05, "*", "")  
          results[[length(results) + 1]] <- data.frame( 
            Pop = "All",  # No grouping, use "All"
            Morph = char,
            shapiro.stat = shapiro.res$statistic,
            shapiro.pvalue = p.value,
            significance = significance
          )
        } else {
          results[[length(results) + 1]] <- data.frame(
            Pop = "All",  
            Morph = char,
            shapiro.stat = NA,
            shapiro.pvalue = NA,
            significance = NA
          )
        }
      }
    }
  } else {
    ## Perform Shapiro-Wilk test grouped by the Pop column
    for (pop in unique(data[[grp]])) {
      pop.data <- subset(data, data[[grp]] == pop)  
      for (char in names(pop.data)) {
        if (is.numeric(pop.data[[char]])) { 
          if (length(unique(pop.data[[char]])) > 1) {  
            shapiro.res <- shapiro.test(pop.data[[char]])  
            p.value <- shapiro.res$p.value  
            significance <- ifelse(p.value < 0.05, "*", "")  
            results[[length(results) + 1]] <- data.frame( 
              Pop = pop,
              Morph = char,
              shapiro.stat = shapiro.res$statistic,
              shapiro.pvalue = p.value,
              significance = significance
            )
          } else {
            results[[length(results) + 1]] <- data.frame(  
              Pop = pop,
              Morph = char,
              shapiro.stat = NA,
              shapiro.pvalue = NA,
              significance = NA
            )
          }
        }
      }
    }
  }

  ## Combine all results into a single data frame and return it
  shapiro.sum <- do.call(rbind, results)
  return(shapiro.sum)
}

