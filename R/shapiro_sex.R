#' @title Shapiro-Wilk Test by Sex
#'
#' @description A function to test normality of data in each sex within each population using Shapiro-Wilk Test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#'
#' shapiro.res <- shapiro_sex(data, sex = "Sex", grp = "Pop")
#'
#' shapiro.res
#'
#' @importFrom stats shapiro.test
#' @return A data frame summarizing the results of Shapiro-Wilk Test for each sex within a population or species.
#' @export

shapiro_sex <- function(data, sex = "Sex", grp = "Pop") {
  ## Check if the specified columns exist in the dataframe
  if (!(sex %in% colnames(data))) {
    stop(paste("Error: Column", sex, "not found in the dataframe."))
  }
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the dataframe."))
  }
  
  ## Extract unique combinations of Pop and Sex
  pop.sex <- unique(data[, c(grp, sex)])
  
  results <- list()

  ## Loop through each unique combination of Pop and Sex
  for (i in 1:nrow(pop.sex)) {
    pop <- pop.sex[[grp]][i]
    sex_value <- pop.sex[[sex]][i]
    subset_data <- subset(data, data[[grp]] == pop & data[[sex]] == sex_value)
    
    ## Perform Shapiro-Wilk test
    shapiro.res <- shapiro.test(subset_data$SVL)
    
    ## Determine significance
    significance <- ifelse(shapiro.res$p.value < 0.05, "*", "")
    
    ## Store the results
    results[[i]] <- data.frame(
      Pop = pop,
      Sex = sex_value,
      shapiro.stat = shapiro.res$statistic,
      shapiro.pvalue = shapiro.res$p.value,
      significance = significance
    )
  }

  ## Combine all results into one data frame
  shapiro.sum <- do.call(rbind, results)
  
  return(shapiro.sum)
}
