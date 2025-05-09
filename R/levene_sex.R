#' @title Levene's Test by Sex
#'
#' @description A function to test homogeneity of variance in each sex within each population using Levene's test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param char String. Column name for morphological character to test.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#'
#' levene.res <- levene_sex(data, sex = "Sex", char = "SVL", grp = "Pop")
#'
#' @importFrom car leveneTest
#' @importFrom stats as.formula
#' @return A data frame summarizing the results of Levene's test for each sex within a population or species.
#' @export

levene_sex <- function(data, sex = "Sex", char = "SVL", grp = "Pop") {
  ## Check if the specified columns exist in the data frame
  if (!(sex %in% colnames(data))) {
    stop(paste("Error: Column", sex, "not found in the data frame."))
  }
  if (!(char %in% colnames(data))) {
    stop(paste("Error: Column", char, "not found in the data frame."))
  }
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the data frame."))
  }

  data$Sex <- as.factor(data$Sex)

  ## Initialize a data frame to store the results
  levene.sum <- data.frame(
    Pop = character(),
    F.value = numeric(),
    p.value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  ## Iterate through each unique population
  for (pop in unique(data[[grp]])) {
    pop.data <- subset(data, data[[grp]] == pop) 
    levene.res <- leveneTest(as.formula(paste(char, "~", sex)), data = pop.data)
    f.val <- levene.res$`F value`[1]
    p.val <- levene.res$`Pr(>F)`[1]  
    significance <- ifelse(p.val < 0.05, "*", "")
    levene.sum <- rbind(levene.sum, data.frame( 
      Pop = pop,
      F.value = f.val,
      p.value = p.val,
      significance = significance
    ))
  }
  
  ## Return the result
  return(levene.sum)  
}
