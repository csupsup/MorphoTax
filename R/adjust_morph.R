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
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' adj.data <- adjust_morph(data, type = "multi_pop")
#' 
#' head(adj.data)
#'
#' @references
#' Thorpe, R. S. (1975). Quantitative handling of characters useful in snake systematics with particular reference to intraspecific 
#' variation in the Ringed Snake \emph{Natrix natrix} (L.). Biological Journal of the Linnaean Society 7: 27–43. https://doi.org/10.1111/j.1095-8312.1975.tb00732.x
#'
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom dplyr bind_rows
#' @return A data frame containing the adjusted morphological data.
#' @export

adjust_morph <- function(data, type = c("multi_pop", "single_pop", "species")) {
    
  ## Option to handle both single population and multiple populations
  type <- match.arg(type)
  
  ## Make an empty list to store adjusted data
  adjusted_data_list <- list()
  
  ## Check if "single population" or "multiple populations" is provided
  if (type == "multi_pop") {
    
    ## Loop over each population
    for (pop in unique(data[, 1])) {
      ## Subset population
      pop_data <- subset(data, data[, 1] == pop) # Pop
      
      ## Calculate the mean of SVL for the population
      SVL_mean <- mean(data$SVL, na.rm = TRUE)

      ## Log-transform SVL
      pop_data$SVL_log <- log10(pop_data$SVL)
      
      ## Loop through all characters (except Pop and SVL)
      numeric_cols <- sapply(pop_data, is.numeric)
      numeric_cols["Pop"] <- FALSE  
      numeric_cols["SVL"] <- FALSE 
      
      for (col in names(pop_data)[numeric_cols]) {
        ## Skip columns with missing values
        if (any(is.na(pop_data[[col]]))) next
        
        ## Calculate the slope (b) using linear regression (log10(X) ~ log10(SVL))
        model <- lm(log10(pop_data[[col]]) ~ log10(pop_data$SVL))
        b <- coef(model)[2]  # slope
        
        ## Apply body size correction: X_adj = log10(X) - b * (log10(SVL) - log10(SVL_mean))
        pop_data[[paste0(col, "_adj")]] <- log10(pop_data[[col]]) - b * (log10(pop_data$SVL) - log10(SVL_mean))
      }
      
      ## Combine adjusted data, SVL_log, and Pop to the list
      adjusted_data <- pop_data[, c("Pop", "SVL_log", grep("_adj", colnames(pop_data), value = TRUE))]
      
      ## Exclude 'SVL_log_adj' if it's included
      adjusted_data <- adjusted_data[, !grepl("SVL_log_adj", colnames(adjusted_data))]
      
      ## Store adjusted data by population
      adjusted_data_list[[as.character(pop)]] <- adjusted_data
    }
    
    ## Convert adjusted data into one dataframe
    final_data <- bind_rows(adjusted_data_list)
    
  } else if (type == "single_pop") {
    ## Obtain data (assuming one population)
    pop_data <- data
    
    ## Calculate the mean of SVL
    SVL_mean <- mean(pop_data$SVL, na.rm = TRUE)
    
    ## Log-transform the SVL
    pop_data$SVL_log <- log10(pop_data$SVL)
    
    ## Loop through all characters
    numeric_cols <- sapply(pop_data, is.numeric)
    numeric_cols["SVL"] <- FALSE  ## Exclude SVL
    
    for (col in names(pop_data)[numeric_cols]) {
      ## Skip columns with missing values
      if (any(is.na(pop_data[[col]]))) next
      
      ## Calculate the slope (b) using linear regression (log10(X) ~ log10(SVL))
      model <- lm(log10(pop_data[[col]]) ~ log10(pop_data$SVL))
      b <- coef(model)[2]  # slope
      
      ## Apply body size correction: X_adj = log10(X) - b * (log10(SVL) - log10(SVL_mean))
      pop_data[[paste0(col, "_adj")]] <- log10(pop_data[[col]]) - b * (log10(pop_data$SVL) - log10(SVL_mean))
    }
    
    ## Combine adjusted data, SVL_log, and Pop to the list
    adjusted_data <- pop_data[, c("SVL_log", grep("_adj", colnames(pop_data), value = TRUE))]
    
    ## Store adjusted data
    adjusted_data_list <- adjusted_data
    
    final_data <- adjusted_data_list
  } else if (type == "species") {
    ## Loop over each population
    for (pop in unique(data[, 1])) {
      ## Subset population
      pop_data <- subset(data, data[, 1] == pop) # Pop
      
      ## Calculate the mean of SVL for the population
      SVL_mean <- mean(pop_data$SVL, na.rm = TRUE)

      ## Log-transform SVL
      pop_data$SVL_log <- log10(pop_data$SVL)
      
      ## Loop through all characters (except Pop and SVL)
      numeric_cols <- sapply(pop_data, is.numeric)
      numeric_cols["Pop"] <- FALSE  
      numeric_cols["SVL"] <- FALSE 
      
      for (col in names(pop_data)[numeric_cols]) {
        ## Skip columns with missing values
        if (any(is.na(pop_data[[col]]))) next
        
        ## Calculate the slope (b) using linear regression (log10(X) ~ log10(SVL))
        model <- lm(log10(pop_data[[col]]) ~ log10(pop_data$SVL))
        b <- coef(model)[2]  # slope
        
        ## Apply body size correction: X_adj = log10(X) - b * (log10(SVL) - log10(SVL_mean))
        pop_data[[paste0(col, "_adj")]] <- log10(pop_data[[col]]) - b * (log10(pop_data$SVL) - log10(SVL_mean))
      }
      
      ## Combine adjusted data, SVL_log, and Pop to the list
      adjusted_data <- pop_data[, c("Pop", "SVL_log", grep("_adj", colnames(pop_data), value = TRUE))]
      
      ## Exclude 'SVL_log_adj' if it's included
      adjusted_data <- adjusted_data[, !grepl("SVL_log_adj", colnames(adjusted_data))]
      
      ## Store adjusted data by population
      adjusted_data_list[[as.character(pop)]] <- adjusted_data
    }
    
    ## Convert adjusted data into one dataframe
    final_data <- bind_rows(adjusted_data_list)
    
  }

  ## Clean column names by removing "_log" and "_adj" from the names
  colnames(final_data) <- gsub("_log$", "", colnames(final_data))
  colnames(final_data) <- gsub("_adj$", "", colnames(final_data))
  
  ## Ensure output is a data frame
  final_data <- as.data.frame(final_data)

  return(final_data)
}
