#' @title Chi-square Test 
#'
#' @description This function performs a Chi-square test, followed by 
#' pairwise comparison for significant results. The test applies to all character columns.
#'
#' @param data A data frame where the first column is population/species and the rest are character with binary data.
#' @param grp String specifying the population column.
#' @param alpha Number. Significance level of p-value.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "binary.samp.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#' 
#' chi.res <- chisq_test(data, grp = "Pop", alpha = 0.05)
#' 
#' head(chi.res)
#' 
#' @importFrom stats chisq.test p.adjust
#' @importFrom utils combn
#' 
#' @return A list containing the contingency table, Chi-square results, and Pairwise comparisons.
#' @export

chisq_test <- function(data, grp = "Pop", alpha = 0.05) {
  ## Convert data to contingency tables
  tbl <- list()
  
  ## Loop through all columns
  for (col in names(data)) {
    if (col != grp && !is.character(data[[col]])) {
      
      ## Create contingency table
      temp_tbl <- table(data[[grp]], data[[col]])
      
      ## Rename column names
      if (all(colnames(temp_tbl) %in% c("0", "1"))) {
        colnames(temp_tbl) <- c("Absent", "Present")
      }
      
      ## Store the table
      tbl[[col]] <- temp_tbl
    }
  }
  
  ## Chi-square results data frame
  chi_results <- data.frame(
    char = character(),
    X.squared = numeric(),
    df = numeric(),
    p.value = numeric(), 
    stringsAsFactors = FALSE
  )
  
  ## Loop through each table and run chi-square test
  for (name in names(tbl)) {
    test <- chisq.test(tbl[[name]])
    pval <- test$p.value
    
    chi_results <- rbind(chi_results, data.frame(
      char = name,
      X.squared = as.numeric(test$statistic),
      df = as.numeric(test$parameter),
      p.value = pval,
      stringsAsFactors = FALSE
    ))
  }
  
  ## Perform pairwise tests for significant variables
  pairwise <- list()
  
  ## Get unique groups
  groups <- unique(data[[grp]])
  
  ## Function to do pairwise comparisons
  pairwise_chisq <- function(var) {
    pairs <- combn(groups, 2, simplify = FALSE)
    results <- data.frame(
      comparison = character(),
      X.squared = numeric(),
      df = numeric(),
      p.value = numeric(),
      p.adj = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (pair in pairs) {
      ## Subset data for two groups
      subdata <- data[data[[grp]] %in% pair, ]
      
      ## Contingency table for two groups
      tab <- table(subdata[[grp]], subdata[[var]])
      
      # Rename columns if 0 and 1
      if (all(colnames(tab) %in% c("0", "1"))) {
        colnames(tab) <- c("Absent", "Present")
      }
      
      test <- chisq.test(tab)
      results <- rbind(results, data.frame(
        comparison = paste(pair[1], pair[2], sep = "-"),
        X.squared = as.numeric(test$statistic),
        df = as.numeric(test$parameter),
        p.value = test$p.value,
        p.adj = NA, 
        stringsAsFactors = FALSE
      ))
    }
    ## Bonferroni correction
    results$p.adj <- p.adjust(results$p.value, method = "bonferroni")
    
    return(results)
  }
  
  ## Run pairwise tests only for significant variables
  sig_vars <- chi_results$char[chi_results$p.value < alpha]
  for (v in sig_vars) {
    pairwise[[v]] <- pairwise_chisq(v)
  }
  
  ## Flatten pairwise list into a single dataframe
  if (length(pairwise) > 0) {
    pairwise_flat <- do.call(rbind, lapply(names(pairwise), function(v) {
      df <- pairwise[[v]]
      df$char <- v
      df
    }))
    rownames(pairwise_flat) <- NULL
    
    ## Reorder columns
    pairwise_flat <- pairwise_flat[, c("char", setdiff(names(pairwise_flat), "char"))]
    
    ## Add significance column
    pairwise_flat$significance <- ifelse(pairwise_flat$p.adj < alpha, "*", "")
  } else {
    pairwise_flat <- data.frame()
  }
  
  ## Return results
  results <- list(
    tbl = tbl,
    chi_results = chi_results,
    pairwise_df = pairwise_flat
  )
  
  return(results)
}
