#' @title Clean P-values
#'
#' @description A function to clean and put together all the p-values from the functions "aov_tukey", "kw_dunn" and mc_chisq_test.
#' All input CSV files should be located in a folder where each file has a prefix "tk_" for Tukey test results, "dn_" for Dunn 
#' and "chi_" for Chi-square tests. Each file contains the results of a test on a certain character, for example, "tk_SVL_results.csv."
#'
#' @param tk.dir A directory that contains the results of the Tukey test.
#' @param dn.dir A directory that stores the results of the Dunn test.
#' @param ch.dir A directory that stores the results of the Chi-square test.
#' @param output.file A directory where the output csv file will be stored.
#' @param asterisk Logical (TRUE or FALSE). If TRUE, it adds an asterisk symbol to a significant p-value.
#'
#' @examples
#' tk_dir <- system.file("extdata", "tk_res", package = "MorphoTax")
#' dn_dir <- system.file("extdata", "dn_res", package = "MorphoTax")
#'
#' clean_pvals(tk.dir = tk_dir, dn.dir = dn_dir, ch.dir = NULL,
#'          output.file = "cleaned_pvals.csv", asterisk = FALSE)
#'
#' @importFrom utils read.csv
#' @importFrom tools file_path_sans_ext
#' @return A data frame containing p-values for all population or species comparisons.
#' @export

clean_pvals <- function(tk.dir = NULL, dn.dir = NULL, ch.dir = NULL,
                        output.file = "cleaned_pvalues.csv", asterisk = TRUE) {

  ## Check if all directories are NULL
  if (is.null(tk.dir) && is.null(dn.dir) && is.null(ch.dir)) {
    stop("Error: At least one of 'tk.dir', 'dn.dir', or 'ch.dir' must be specified.")
  }

  ## Check if specified directories exist
  for (dir in list(tk.dir, dn.dir, ch.dir)) {
    if (!is.null(dir) && !dir.exists(dir)) {
      stop(paste("Error: The directory", dir, "does not exist."))
    }
  }

  ## Clean directories for processed files
  for (directory_path in list(tk.dir, dn.dir, ch.dir)) {
    if (!is.null(directory_path)) {
      files <- list.files(path = directory_path, pattern = "_processed.csv$", full.names = TRUE)
      if (length(files) > 0) file.remove(files)
    }
  }

  ## Function to normalize comparison labels
  norm_comparison <- function(comparison) {
    parts <- strsplit(gsub(" ", "", comparison), "-")
    sorted.parts <- sapply(parts, function(x) paste(sort(x), collapse = "-"))
    return(sorted.parts)
  }

  ## Function to process and write CSVs
  write_files <- function(directory) {
    file.list <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
    lapply(file.list, function(file) {
      processed_file <- file.path(directory, paste0(tools::file_path_sans_ext(basename(file)), "_processed.csv"))
      data <- read.csv(file)
      if ("Comparison" %in% colnames(data)) {
        data$Normalized_Comparison <- norm_comparison(data$Comparison)
      }
      write.csv(data, processed_file, row.names = FALSE)
      cat("Processed file written (or overwritten):", processed_file, "\n")
    })
    invisible(NULL)
  }

  ## Process provided directories
  if (!is.null(tk.dir)) write_files(tk.dir)
  if (!is.null(dn.dir)) write_files(dn.dir)
  if (!is.null(ch.dir)) write_files(ch.dir)

  ## Initialize data frames
  tk.pval <- dn.pval <- ch.pval <- comparisons <- NULL

  # Tukey
  if (!is.null(tk.dir)) {
    tk.tests <- list.files(tk.dir, pattern = '_processed.csv$', full.names = TRUE)
    if (length(tk.tests) > 0) {
      tk.pval <- as.data.frame(lapply(tk.tests, function(x) read.csv(x)$p.adj))
      names(tk.pval) <- tools::file_path_sans_ext(basename(tk.tests))
      colnames(tk.pval) <- gsub("^tk_(.*?)_results.*", "\\1", colnames(tk.pval))

      df <- read.csv(tk.tests[1])
      if ("Normalized_Comparison" %in% colnames(df)) {
        comparisons <- df$Normalized_Comparison
      } else if ("Comparison" %in% colnames(df)) {
        comparisons <- norm_comparison(df$Comparison)
      }
    }
  }

  # Dunn
  if (!is.null(dn.dir)) {
    dn.tests <- list.files(dn.dir, pattern = '_processed.csv$', full.names = TRUE)
    if (length(dn.tests) > 0) {
      dn.pval <- as.data.frame(lapply(dn.tests, function(x) read.csv(x)$P.adj))
      names(dn.pval) <- tools::file_path_sans_ext(basename(dn.tests))
      colnames(dn.pval) <- gsub("^dn_(.*?)_results.*", "\\1", colnames(dn.pval))
    }
  }

  # Chi-square
  if (!is.null(ch.dir)) {
    ch.tests <- list.files(ch.dir, pattern = '^chi_.*?_results_processed\\.csv$', full.names = TRUE)
    if (length(ch.tests) > 0) {
      ch.pval <- as.data.frame(lapply(ch.tests, function(x) read.csv(x)$p.adj))
      names(ch.pval) <- tools::file_path_sans_ext(basename(ch.tests))
      colnames(ch.pval) <- gsub("^chi_(.*?)_results.*", "\\1", colnames(ch.pval))
    }
  }

  ## Combine all p-values
  if (is.null(comparisons)) {
    stop("Could not determine comparisons. Ensure at least one Tukey file with comparisons is present.")
  }

  phoc.pval <- data.frame(comparisons)
  if (!is.null(tk.pval)) phoc.pval <- cbind(phoc.pval, tk.pval)
  if (!is.null(dn.pval)) phoc.pval <- cbind(phoc.pval, dn.pval)
  if (!is.null(ch.pval)) phoc.pval <- cbind(phoc.pval, ch.pval)

  ## Round and add asterisks
  phoc.pval[] <- lapply(phoc.pval, function(x) {
    if (is.numeric(x)) {
      x_rounded <- round(x, 3)
      if (asterisk) {
        ifelse(x_rounded < 0.05, paste0(format(x_rounded, nsmall = 3), "*"), format(x_rounded, nsmall = 3))
      } else {
        format(x_rounded, nsmall = 3)
      }
    } else {
      x
    }
  })

  ## Convert all columns to numeric except the first one
  phoc.pval[-1] <- lapply(phoc.pval[-1], as.numeric)

  ## Write output
  if (!is.null(output.file)) {
    write.csv(phoc.pval, output.file, row.names = FALSE)
    cat("Processed p-values written to:", output.file, "\n")
  } else {
    cat("Output file not specified. P-values not written.\n")
  }

  return(phoc.pval)
}
