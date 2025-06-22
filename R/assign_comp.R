#' @title Assign Comparisons
#'
#' @description A function for assigning groups (e.g., species) to each comparison in a dataframe of p-values, 
#' allowing the comparisons to be sorted by groups using the function 'plot_pval_bubble.'
#' 
#' @param data A data frame with p-values of comparisons.
#' @param comp String specifying the column name of comaparisons
#' @param rules A list containing the group assignment rules.
#' 
#' @examples
#' data <- read.csv(system.file("extdata", "cleaned.pvals.csv", package = "MorphoTax"))
#' sp.grp <- list(
#'            marmorata = c("Luzon"),
#'            philippina = c("Mindanao"),
#'            mallarii = c("Palawan"))
#' 
#' data <- assign_comp(data, comp = "comparisons", rules = sp.grp)
#'
#' @return A dataframe with group assignment column.
#' @export

assign_comp <- function(data, comp = "comparisons", rules) {
  ## Internal function to detect species based on rules
  detect_species <- function(comparison, rules) {
    matched_species <- c()
    for (species in names(rules)) {
      for (pattern in rules[[species]]) {
        if (grepl(pattern, comparison)) {
          matched_species <- unique(c(matched_species, species))
          break
        }
      }
    }
    if (length(matched_species) > 0) {
      return(paste(sort(matched_species), collapse = "+"))
    } else {
      return(NA_character_)
    }
  }

  ## Check that the comparison column exists
  if (!comp %in% names(data)) {
    stop(paste("Column", comp, "not found in data."))
  }

  ## Create species vector
  species_vec <- sapply(data[[comp]], detect_species, rules = rules)

  ## Insert species column after the comparison column
  comp_index <- match(comp, names(data))
  data <- cbind(
    data[1:comp_index],
    species = species_vec,
    data[(comp_index + 1):ncol(data)]
  )

  rownames(data) <- NULL

  return(data)
}
