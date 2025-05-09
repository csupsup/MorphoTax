#' @title Levene's Test
#'
#' @description A function to test homogeneity of variance in all morphological data using Levene's test.
#' 
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#'
#' levene.res <- levene_all(data, grp = "Pop")
#'
#' @importFrom car leveneTest
#' @return A data frame summarizing the results of Levene's test.
#' @export

levene_all <- function(data, grp = "Pop") {
  ## Check if the input 'data' is a data frame
  if (!is.data.frame(data)) {
    stop("Error: The input data is not a data frame.")
  }
  
  ## Ensure the group column is a factor
  data[[grp]] <- as.factor(data[[grp]])

  ## Initialize an empty data frame to store results
  levene.sum <- data.frame(
    Morph = character(),
    F.value = numeric(),
    p.value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )

  ## Loop through each column in the data
  for (char in colnames(data)) {
    # Check if the column is numeric
    if (is.numeric(data[[char]])) {
      # Perform Levene's test
      levene.res <- leveneTest(data[[char]] ~ data[[grp]])
      
      ## Extract the F-value and p-value from the result
      levene.stat <- levene.res$`F value`[1]
      p.val <- levene.res$`Pr(>F)`[1]
      
      ## Determine the significance
      significance <- ifelse(p.val < 0.05, "*", "")
      
      ## Add the results to the summary data frame
      levene.sum <- rbind(levene.sum, data.frame(
        Morph = char,
        F.value = levene.stat,
        p.value = p.val,
        significance = significance
      ))
    }
  }

  ## Return the summary data frame
  return(levene.sum)
}
