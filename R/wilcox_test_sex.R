#' @title Wilcoxon Test by Sex
#'
#' @description A function to test differences in characters between sexes within each population using Wilcoxon Test.
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param sex String. Column name for sex information, specifying male and female.
#' @param char String. Column name for morphological character to test.
#' @param grp String. Column name for population or species.
#' @param all Logical. If TRUE, the Wilcoxon Test will be performed to all character or numeric columns.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#'
#' wilcox.res <- wilcox_test_sex(data, sex = "Sex", char = "SVL", grp = "Pop", all = FALSE)
#'
#' head(wilcox.res)
#'
#' @importFrom stats wilcox.test
#' @return A data frame containing the results of Wilcoxon test.
#' @export

wilcox_test_sex <- function(data, sex = "Sex", char = "SVL", grp = "Pop", all = FALSE) {
  ## Check if the specified columns exist in the data frame
  if (!(sex %in% colnames(data))) {
    stop(paste("Error: Column", sex, "not found in the data frame."))
  }
  if (!(grp %in% colnames(data))) {
    stop(paste("Error: Column", grp, "not found in the data frame."))
  }

  ## Determine which columns to analyze
  if (all) {
    num.cols <- sapply(data, is.numeric)
    char.cols <- setdiff(names(data)[num.cols], c(grp, sex))
    if (length(char.cols) == 0) {
      stop("No numeric columns to analyze.")
    }
  } else {
    if (!(char %in% colnames(data))) {
      stop(paste("Error: Column", char, "not found in the data frame."))
    }
    char.cols <- char
  }

  ## Initialize result list
  result.list <- list()

  ## Loop through each character column
  for (char.col in char.cols) {
    test.sum <- data.frame(
      Char = character(),
      Pop = character(),
      W.stat = numeric(),
      p.value = numeric(),
      significance = character(),
      stringsAsFactors = FALSE
    )

    for (pop in unique(data[[grp]])) {
      pop.data <- subset(data, data[[grp]] == pop)

      ## Ensure both sex groups are present
      if (length(unique(pop.data[[sex]])) < 2) {
        warning(paste("Skipping population", pop, "- less than two groups in", sex))
        next
      }

      test.res <- wilcox.test(as.formula(paste(char.col, "~", sex)), data = pop.data, exact = FALSE)

      test.sum <- rbind(test.sum, data.frame(
        char = char.col,
        pop = pop,
        W.stat = test.res$statistic,
        p.value = test.res$p.value,
        significance = ifelse(test.res$p.value < 0.05, "*", "")
      ))
    }

    result.list[[char.col]] <- test.sum
  }

  ## Combine results into a single data frame
  final.result <- if (length(result.list) == 1) {
    result.list[[1]]
  } else {
    do.call(rbind, result.list)
  }

  row.names(final.result) <- NULL
  return(final.result)
}
