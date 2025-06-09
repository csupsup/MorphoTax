#' @title Levene's Test by Sex
#'
#' @description A function to test homogeneity of variance in each sex within each population using Levene's test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#'
#' levene.res <- levene_sex(data, sex = "Sex", grp = "Pop")
#'
#' levene.res
#' 
#' @importFrom car leveneTest
#' @importFrom stats as.formula
#' @return A data frame summarizing the results of Levene's test for each sex within a population or species.
#' @export

levene_sex <- function(data, sex = "Sex", grp = "Pop") {
  ## Check if the specified columns exist in the data frame
  if (!(sex %in% colnames(data))) {
    stop(paste("Error: Column", sex, "not found in the data frame."))
  }
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the data frame."))
  }
  
  data[[sex]] <- as.factor(data[[sex]])
  
  ## Identify numeric columns
  numeric_cols <- sapply(data, is.numeric)
  numeric_names <- names(data)[numeric_cols]
  
  ## Initialize a data frame to store the results
  levene.sum <- data.frame(
    Pop = character(),
    Char = character(),
    F.value = numeric(),
    p.value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  ## Iterate through each unique population
  for (pop in unique(data[[grp]])) {
    pop.data <- subset(data, data[[grp]] == pop)
    
    ## Loop through each numeric column
    for (var in numeric_names) {
      levene.res <- leveneTest(as.formula(paste(var, "~", sex)), data = pop.data)
      f.val <- levene.res$`F value`[1]
      p.val <- levene.res$`Pr(>F)`[1]
      significance <- ifelse(p.val < 0.05, "*", "")
      
      levene.sum <- rbind(levene.sum, data.frame(
        Pop = pop,
        Char = var,
        F.value = f.val,
        p.value = p.val,
        significance = significance
      ))
    }
  }
  
  ## Return the result
  rownames(levene.sum) <- NULL
  
  return(levene.sum)
}