#' @title Filter Data by Sex
#'
#' @description A function to filter data by excluding populations represented by only one sex.
#' Additionally, it filters sex data within each population with fewer than three samples.
#' 
#' @param data A dataframe with population and sex information in the first two columns, followed by morphology data.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' 
#' filtered.data <- filter_sex(data) 
#'
#' head(filtered.data)
#'
#' @return A data frame containing filtered data.
#' @export

filter_sex <- function(data) {
  ## Ensure "Sex" column is a factor
  data$Sex <- as.factor(data$Sex)
  
  ## Display initial sex counts
  sex.counts <- table(data$Pop, data$Sex)
  cat("Initial sex counts:\n")
  print(sex.counts)
  
  ## Remove populations with only one sex represented
  sex.occ <- table(data$Pop, data$Sex)
  sex.comp <- rownames(sex.occ)[rowSums(sex.occ > 0) > 1]
  removed_one_sex <- setdiff(rownames(sex.occ), sex.comp)
  data <- data[data$Pop %in% sex.comp, ]
  
  ## Display removed populations
  if (length(removed_one_sex) > 0) {
    cat("\nPopulations removed due to having only one sex represented:\n")
    print(removed_one_sex)
  } else {
    cat("\nNo populations removed for having only one sex represented.\n")
  }
  
  ## Display filtered sex counts after removing populations with only one sex
  sex.counts <- table(data$Pop, data$Sex)
  cat("\nFiltered sex counts (after removing populations with only one sex represented):\n")
  print(sex.counts)
  
  ## Filter populations where each sex has >4 samples
  sex.counts <- table(data$Pop, data$Sex)
  check.pops <- rownames(sex.counts)[apply(sex.counts, 1, function(x) all(x >= 4))]
  removed_few_samples <- setdiff(rownames(sex.counts), check.pops)
  data <- data[data$Pop %in% check.pops, ]
  
  ## Display removed populations
  if (length(removed_few_samples) > 0) {
    cat("\nPopulations removed due to having fewer than 4 samples for each sex:\n")
    print(removed_few_samples)
  } else {
    cat("\nNo populations removed for having fewer than 4 samples for each sex.\n")
  }
  
  ## Display filtered sex counts after the second filtering step
  sex.counts <- table(data$Pop, data$Sex)
  cat("\nFiltered sex counts (after filtering populations with less than 4 samples per sex):\n")
  print(sex.counts)
  
  ## Return the filtered dataframe
  return(data)
}

