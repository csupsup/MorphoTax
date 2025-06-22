#' @title Chi-square Test with Monte Carlo
#'
#' @description This function performs a Monte Carlo Chi-square test, followed by 
#' pairwise Fisher’s exact tests for significant results. The test applies to all character columns.
#'
#' @param data A data frame where the first column is population/species and the rest are character with binary data.
#' @param grp String specifying the population column.
#' @param alpha Number. Significance level of p-value.
#' @param write.res Logical. If TRUE, writes the summary results of every character to a csv file.
#' @param dir String specifying the output directory.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "binary.samp.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#' 
#' mc.chi.res <- mc_chisq_test(data, grp = "Pop", alpha = 0.05, write.res = FALSE, dir = NULL)
#' 
#' head(mc.chi.res)
#' 
#' @importFrom rcompanion pairwiseNominalIndependence
#' @importFrom dplyr mutate rename select filter
#' @importFrom stats chisq.test 
#' 
#' @return A list containing the contingency table, Monte Carlo Chi-square results, and Pairwise Fisher’s exact tests.
#' @export

mc_chisq_test <- function(data, grp = "Pop", alpha = 0.05, write.res = FALSE, dir = NULL) {
  if (!requireNamespace("rcompanion", quietly = TRUE)) {
    stop("Please install the 'rcompanion' package: install.packages('rcompanion')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package: install.packages('dplyr')")
  }

  ## Initialize storage
  tbl <- list()
  pairwise_results <- list()

  ## Create contingency tables
  for (col in names(data)) {
    if (col != grp && !is.character(data[[col]])) {
      temp_tbl <- table(data[[col]], data[[grp]])

      ## Rename rows if binary 0/1
      if (all(rownames(temp_tbl) %in% c("0", "1"))) {
        rownames(temp_tbl) <- c("Absent", "Present")
      }

      ## Transpose
      temp_tbl <- t(temp_tbl)

      tbl[[col]] <- temp_tbl
    }
  }

  ## Chi-squared results data frame
  chi_results <- data.frame(
    char = character(),
    X.squared = numeric(),
    df = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )

  ## Run tests
  for (name in names(tbl)) {
    test <- tryCatch({
      chisq.test(tbl[[name]], simulate.p.value = TRUE, B = 1e5)
    }, warning = function(w) {
      chisq.test(tbl[[name]])
    }, error = function(e) {
      list(statistic = NA, parameter = NA, p.value = NA)
    })

    pval <- as.numeric(test$p.value)
    df_val <- (nrow(tbl[[name]]) - 1) * (ncol(tbl[[name]]) - 1)

    chi_results <- rbind(chi_results, data.frame(
      char = name,
      X.squared = ifelse(is.null(test$statistic), NA, round(as.numeric(test$statistic), 5)),
      df = df_val,
      p.value = ifelse(is.na(pval), NA, round(pval, 5)),
      stringsAsFactors = FALSE
    ))

    ## If significant, run pairwise test
    if (!is.na(pval) && pval < alpha) {
      pw <- tryCatch({
        pairwiseNominalIndependence(tbl[[name]],
                                    fisher = TRUE,
                                    gtest = FALSE,
                                    chisq = FALSE,
                                    method = "holm")
      }, error = function(e) {
        warning(paste("Pairwise test failed for", name, ":", e$message))
        NULL
      })

      if (!is.null(pw)) {
        pairwise_results[[name]] <- pw
      }
    }
  }

  ## Convert pairwise results to a data frame
  pairwise_summary_df <- function(pairwise_list, alpha_val) {
    all_results <- lapply(names(pairwise_list), function(varname) {
      df <- pairwise_list[[varname]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df %>%
        mutate(char = varname) %>%
        rename(Comparison = Comparison,
               p.value = p.Fisher,
               p.adj = p.adj.Fisher) %>%
        mutate(significance = ifelse(p.adj < alpha_val, "*", "")) %>%
        select(char, Comparison, p.value, p.adj, significance)
    })
    combined_df <- do.call(rbind, all_results)
    rownames(combined_df) <- NULL
    return(combined_df)
  }

  ## Generate the final pairwise df
  pairwise_df <- pairwise_summary_df(pairwise_results, alpha)
  if (!is.null(pairwise_df)) {
    pairwise_df$Comparison <- gsub(" : ", "-", pairwise_df$Comparison)
  }

  ## Optionally write results to CSV
  if (write.res && !is.null(pairwise_df)) {
    ## Create output directory only if specified
    if (!is.null(dir) && !dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }

    unique_chars <- unique(pairwise_df$char)
    for (char in unique_chars) {
      char_df <- filter(pairwise_df, char == !!char)
      names(char_df)[names(char_df) == "char"] <- "Morph"

      if (!is.null(dir)) {
        fname <- file.path(dir, paste0("chi_", char, "_results.csv"))
      } else {
        fname <- paste0("chi_", char, "_results.csv")
      }

      write.csv(char_df, file = fname, row.names = FALSE)
    }
  }

  return(list(
    tbl = tbl,
    chi_results = chi_results,
    pairwise = pairwise_results,
    pairwise_df = pairwise_df
  ))
}


## Declare
utils::globalVariables(c("p.Fisher", "p.adj.Fisher", "char", "Comparison", "p.value", "p.adj", "significance"))