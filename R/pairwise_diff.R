#' @title Pairwise Relative Difference
#' 
#' @description A function to calculate pairwise relative differences between population means of a specific morphological character. 
#' Calculated as:
#' 
#' rd = abs(pop1Mean - pop2Mean) / ((pop1Mean + pop2Mean) / 2) * 100
#'
#' @param data A data frame with population or species label in the first column, followed by SVL and other data.
#' @param grp String indicating the column name of the population.
#' @param char String specifying the column name of the character.
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' 
#' pw.diff <- pairwise_diff(data, grp = "Pop", char = "SVL")
#' 
#' head(pw.diff)
#' 
#' @importFrom dplyr group_by summarise rename mutate left_join case_when
#' 
#' @return A data frame with all pairwise comparisons between populations,
#' including relative differences and which population had the higher mean.
#'
#' @export

pairwise_diff <- function(data, grp = "Pop", char = "TL") {

  ## Calculate mean character per population
  mean_char <- data %>%
    group_by(.data[[grp]]) %>%
    summarise(mean_val = mean(.data[[char]], na.rm = TRUE), .groups = "drop")
  
  ## Create all pairwise combinations excluding self comparisons
  pairs <- expand.grid(mean_char[[grp]], mean_char[[grp]], stringsAsFactors = FALSE)
  colnames(pairs) <- c("pop1", "pop2")
  pairs <- pairs[pairs$pop1 != pairs$pop2, ]
  
  ## Combine means to pairs
  pairs <- pairs %>%
    left_join(mean_char, by = c("pop1" = grp)) %>%
    rename(pop1Mean = mean_val) %>%
    left_join(mean_char, by = c("pop2" = grp)) %>%
    rename(pop2Mean = mean_val) %>%
    mutate(
      relDiff = abs(pop1Mean - pop2Mean) / ((pop1Mean + pop2Mean) / 2),
      relDiffPct = relDiff * 100,
      
      ## Adjust comment name
      !!paste0("longer_", char) := case_when(
        pop1Mean > pop2Mean ~ pop1,
        pop1Mean < pop2Mean ~ pop2,
        TRUE ~ "Equal"
      )
    )
  
  return(pairs)
}

## Declare
utils::globalVariables(c(":=", "pop1Mean", "pop2Mean", "mean_val", "relDiff"))
