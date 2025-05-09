#' @title T-test by Sex
#'
#' @description A function to test differences in characters between sexes with each population using independent T-test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param char String. Column name for morphological character to test.
#' @param grp String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "morphR"))
#'
#' t.test.res <- t_test_sex(data, sex = "Sex", char = "SVL", grp = "Pop")
#'
#' head(t.test.res)
#'
#' @importFrom stats t.test
#' @return A data frame containing the results of T-test.
#' @export

t_test_sex <- function(data, sex = "Sex", char = "SVL", grp = "Pop") {
  # Check if the specified columns exist in the data frame
  if (!(sex %in% colnames(data))) {
    stop(paste("Error: Column", sex, "not found in the data frame."))
  }
  if (!(char %in% colnames(data))) {
    stop(paste("Error: Column", char, "not found in the data frame."))
  }
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the data frame."))
  }

  # Initialize a data frame to store the results
  t.test.sum <- data.frame(
    Pop = character(),
    T.stat = numeric(),
    F.value = numeric(),
    p.value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  ## Iterate through each unique population
  for (pop in unique(data[[grp]])) {
    pop.data <- subset(data, data[[grp]] == pop) 
    t.test.res <- t.test(as.formula(paste(char, "~", sex)), data = pop.data)  
    t.stat <- t.test.res$statistic[1]  
    f.val <- t.test.res$parameter[1]  
    p.val <- t.test.res$p.value  
    significance <- ifelse(p.val < 0.05, "*", "")  
    t.test.sum <- rbind(t.test.sum, data.frame( 
      Pop = pop,
      T.stat = t.stat,
      F.value = f.val,
      p.value = p.val,
      significance = significance
    ))
  }
  
  return(t.test.sum)  
}
