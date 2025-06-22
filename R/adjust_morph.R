#' @title Remove Allometric Effects
#'
#' @description A function to remove the allometric effects of body size in morphological data.
#' Allometric size correction can be applied to both single and multi-population data, as well as species data.
#' The \emph{"single_pop"} option considers all data from one species; therefore, an overall mean is utilized to calculate a single slope for the entire dataset.
#' The \emph{"multi_pop"} option considers all data from one species, while the slope is calculated for each population using an overall mean.
#' The \emph{"species"} option considers each unique group as a single species; therefore, the slope is calculated for each using their respective mean. 
#'
#' @param data A data frame with population or species label in the first column, followed by SVL and other data.
#' @param type String. Data type:  "single_pop", "multi_pop" and "species".
#' @param char String. A character as size reference (e.g., SVL).
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' adj.data <- adjust_morph(data, type = "multi_pop", char = "SVL")
#' 
#' head(adj.data)
#'
#' @references
#' Thorpe, R. S. (1975). Quantitative handling of characters useful in snake systematics with particular reference to intraspecific 
#' variation in the Ringed Snake \emph{Natrix natrix} (L.). Biological Journal of the Linnaean Society 7: 27–43. https://doi.org/10.1111/j.1095-8312.1975.tb00732.x
#'
#' @importFrom stats lm
#' @importFrom stats coef
#' 
#' @return A data frame containing the adjusted morphological data.
#' @export

adjust_morph <- function(data, type = c("multi_pop", "single_pop", "species"), char = "SVL") {
  type <- match.arg(type)

  ## Ensure "Pop" column exists
  if (!"Pop" %in% names(data)) {
    names(data)[1] <- "Pop"
  }

  ## Track original row order
  data$.__row_id <- seq_len(nrow(data))

  ## Global mean of the reference character
  global_char_mean <- mean(data[[char]], na.rm = TRUE)

  ## Adjust one group (species or population)
  adjust_group <- function(group_data, use_global_mean = FALSE) {
    char_log <- log10(group_data[[char]])
    group_data[[paste0(char, "_log")]] <- char_log

    ## Identify which mean to use
    char_mean <- if (use_global_mean) global_char_mean else mean(group_data[[char]], na.rm = TRUE)

    ## Find numeric traits to adjust
    numeric_cols <- sapply(group_data, is.numeric)
    exclude_cols <- c(".__row_id", "Pop", char)
    trait_cols <- setdiff(names(group_data)[numeric_cols], exclude_cols)

    for (trait in trait_cols) {
      if (anyNA(group_data[[trait]])) next
      model <- lm(log10(group_data[[trait]]) ~ char_log)
      b <- coef(model)[2]
      group_data[[paste0(trait, "_adj")]] <- log10(group_data[[trait]]) - b * (char_log - log10(char_mean))
    }

    ## Return relevant columns
    keep_cols <- c(".__row_id", "Pop", paste0(char, "_log"), grep("_adj$", names(group_data), value = TRUE))
    return(group_data[, keep_cols, drop = FALSE])
  }

  ## Apply based on type
  if (type == "multi_pop") {
    ## Use global mean but group-specific slope
    final_data <- do.call(rbind, lapply(split(data, data$Pop), adjust_group, use_global_mean = TRUE))

  } else if (type == "species") {
    ## Use species-specific mean and slope
    final_data <- do.call(rbind, lapply(split(data, data$Pop), adjust_group, use_global_mean = FALSE))

  } else if (type == "single_pop") {
    if (!"Pop" %in% names(data)) data$Pop <- "Population_1"
    final_data <- adjust_group(data, use_global_mean = TRUE)
  }

  ## Rename columns (preserve Pop name)
  new_names <- names(final_data)
  new_names <- ifelse(grepl(paste0("^", char, "_log$"), new_names),
                      char, gsub("_adj$", "", new_names))
  names(final_data) <- new_names

  ## Restore original order
  final_data <- final_data[order(final_data$.__row_id), ]
  final_data$.__row_id <- NULL
  final_data[[paste0(char, "_log")]] <- NULL

  return(as.data.frame(final_data))
}
