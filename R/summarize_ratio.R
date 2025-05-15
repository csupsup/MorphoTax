#' @title Summarize Character Ratio
#'
#' @description A function to calculate and summarize the ratio of two morphological data. The summary includes the number of samples with longer length
#' (ratio >1.0) and shorter length (ratio <1.0), percentages of samples relative to the total, range (min–max), mean, and the calculated 
#' relative difference (e.g., |ED - SVL|/SV) * 100). 
#' 
#' @param data A data frame with population or species label and sex information in the first two columns, followed by morpholigcal data.
#' @param char.as.num String. Column name for character representing the numerator.
#' @param char.as.den String. Column name for character representing the denominator.
#' @param grp  String. Column name for population or species.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#'
#' ratio.sum <- summarize_ratio(data, char.as.num = "TD", char.as.den = "SVL", grp = "Pop")
#'
#' head(ratio.sum)
#'
#' @importFrom stats sd
#' @return A data frame containing summary of ratio.
#' @export

summarize_ratio <- function(data, char.as.num, char.as.den, grp) {
    ## Check if the columns are present in the data frame
    missing_cols <- setdiff(c(char.as.num, char.as.den, grp), colnames(data))
    if (length(missing_cols) > 0) {
        stop(paste("The following columns are missing from the data frame:", paste(missing_cols, collapse = ", ")))
    }

    ## Create the ratio column
    data$ratio <- data[[char.as.num]] / data[[char.as.den]]

    ## Create an empty data frame to store results
    ratio.count <- data.frame(
        Pop = character(0),
        N = numeric(0),
        LongCount = numeric(0),
        ShortCount = numeric(0),
        LongPct = numeric(0),
        ShortPct = numeric(0),
        LongRange = character(0),
        ShortRange = character(0),
        meanLongRatio = numeric(0),
        sdLongRatio = numeric(0),
        meanShortRatio = numeric(0),
        sdShortRatio = numeric(0),
        ShortRdiff = character(0), 
        LongRdiff = character(0) 
    )

    ## Unique populations
    populations <- unique(data[[grp]])

    for (pop in populations) {
        pop_data <- subset(data, data[[grp]] == pop)
        total_samples <- nrow(pop_data)

        above_1 <- sum(pop_data$ratio > 1)
        below_1 <- sum(pop_data$ratio < 1)

        above_1_percent <- (above_1 / total_samples) * 100
        below_1_percent <- (below_1 / total_samples) * 100

        mean_above_1_ratio <- mean(pop_data$ratio[pop_data$ratio > 1], na.rm = TRUE)
        sd_above_1_ratio <- sd(pop_data$ratio[pop_data$ratio > 1], na.rm = TRUE)

        mean_below_1_ratio <- mean(pop_data$ratio[pop_data$ratio < 1], na.rm = TRUE)
        sd_below_1_ratio <- sd(pop_data$ratio[pop_data$ratio < 1], na.rm = TRUE)

        below_1_percent_comment <- ""
        above_1_percent_comment <- "" 

        if (total_samples == 1) {
            above_1_range <- as.character(pop_data$ratio[pop_data$ratio > 1])
            below_1_range <- as.character(pop_data$ratio[pop_data$ratio < 1])
        } else {
            if (above_1 > 0) {
                min_above <- round(min(pop_data$ratio[pop_data$ratio > 1]), 2)
                max_above <- round(max(pop_data$ratio[pop_data$ratio > 1]), 2)
                above_1_range <- paste0(min_above, " to ", max_above)

                ## Add percentage-based comment for longer ratios
                min_pct_longer <- round((min_above - 1) * 100)
                max_pct_longer <- round((max_above - 1) * 100)
                above_1_percent_comment <- paste0(char.as.num, " is ", min_pct_longer, "% to ", max_pct_longer, "% longer than ", char.as.den)
            } else {
                above_1_range <- "NA"
                above_1_percent_comment <- "NA"
            }

            if (below_1 > 0) {
                min_below <- round(min(pop_data$ratio[pop_data$ratio < 1]), 2)
                max_below <- round(max(pop_data$ratio[pop_data$ratio < 1]), 2)
                below_1_range <- paste0(max_below, " to ", min_below)

                ## Add percent-based comment for shorter ratios
                min_pct <- round((1 - min_below) * 100)
                max_pct <- round((1 - max_below) * 100)
                below_1_percent_comment <- paste0(char.as.num, " is ", max_pct, "% to ", min_pct, "% shorter than ", char.as.den)
            } else {
                below_1_range <- "NA"
                below_1_percent_comment <- "NA"
            }
        }

        ratio.count <- rbind(ratio.count, data.frame(
            Pop = pop,
            N = total_samples,
            LongCount = above_1,
            ShortCount = below_1,
            LongPct = above_1_percent,
            ShortPct = below_1_percent,
            LongRange = above_1_range,
            ShortRange = below_1_range,
            meanLongRatio = mean_above_1_ratio,
            sdLongRatio = sd_above_1_ratio,
            meanShortRatio = mean_below_1_ratio,
            sdShortRatio = sd_below_1_ratio,
            ShortRdiff = below_1_percent_comment,  
            LongRdiff = above_1_percent_comment
        ))
    }

    return(ratio.count)
}
