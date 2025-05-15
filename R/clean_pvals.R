#' @title Clean P-values
#'
#' @description A function to clean and put together all the p-values from the functions "aov_tukey" and "kw_dunn".
#' All input CSV files should be located in a folder where each file has a prefix "tk_" for Tukey test results and "dn_" for Dunn test.
#' Each file contains the results of a test on a certain character, for example, "tk_SVL_results.csv."
#'
#' @param tk.dir A directory that contains the results of the Tukey test.
#' @param dn.dir A directory that stores the results of the Dunn test.
#' @param output.file A directory where the output csv file will be stored.
#' @param asterisk Logical (TRUE or FALSE). If TRUE, it adds an asterisk symbol to a significant p-value.
#'
#' @examples
#' tk_dir <- system.file("extdata", "tk_res", package = "MorphoTax")
#' dn_dir <- system.file("extdata", "dn_res", package = "MorphoTax")
#'
#' clean_pvals(tk.dir = tk_dir, dn.dir = dn_dir, 
#'          output.file = "cleaned_pvals.csv", asterisk = FALSE)
#'
#' @importFrom utils read.csv
#' @importFrom tools file_path_sans_ext
#' @return A data frame containing p-values for all population or species comparisons.
#' @export

clean_pvals <- function(tk.dir = "tk_path", dn.dir = "dn_path", 
                       output.file = "cleaned_pvalues.csv", asterisk = TRUE) {

  ## Check if required input directories are provided
  if (missing(tk.dir) || missing(dn.dir)) {
    stop("Error: Both 'tk.dir' and 'dn.dir' must be specified.")
  }
  
  # Check if the directories exist
  if (!dir.exists(tk.dir)) {
    stop(paste("Error: The directory", tk.dir, "does not exist."))
  }
  if (!dir.exists(dn.dir)) {
    stop(paste("Error: The directory", dn.dir, "does not exist."))
  }

  ## Clean directories for processed files
  directories <- c(tk.dir, dn.dir)
  
  ## Iterate over each directory
  for (directory_path in directories) {
    # List all files in the current directory that match the pattern
    files <- list.files(path = directory_path, pattern = "_processed.csv$", full.names = TRUE)

    ## Check if any files match the pattern
    if (length(files) > 0) {
      ## Remove the matching files
      file.remove(files)
    }
  }

  ## Function to check common comparison labels
  norm_comparison <- function(comparison) {
    parts <- strsplit(gsub(" ", "", comparison), "-") 
    sorted.parts <- sapply(parts, function(x) paste(sort(x), collapse = "-"))
    return(sorted.parts)
  }

  ## Function to write processed files (from Tukey or Dunn test results)
  write_files <- function(directory) {
    file.list <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
    
    ## Process each file
    lapply(file.list, function(file) {
      ## Define processed file name
      processed_file <- file.path(directory, paste0(file_path_sans_ext(basename(file)), "_processed.csv"))
      
      ## If processed file already exists, overwrite it
      data <- read.csv(file)
      data$Normalized_Comparison <- norm_comparison(data$Comparison)
      
      ## Always overwrite processed file if it already exists
      write.csv(data, processed_file, row.names = FALSE)
      cat("Processed file written (or overwritten):", processed_file, "\n")
    })
    invisible(NULL)
  }

  ## Process and write files from both directories
  write_files(tk.dir)
  write_files(dn.dir)

  ## Merge p-values from Tukey and Dunn tests
  
  ## Get p-values from Tukey test
  tk.tests <- list.files(tk.dir, pattern = '_processed.csv', full.names = TRUE)
  tk.pval <- as.data.frame(lapply(tk.tests, function(x) read.csv(x)$p.adj))
  names(tk.pval) <- file_path_sans_ext(basename(tk.tests))
  
  ## Fix column names for Tukey p-values
  colnames(tk.pval) <- gsub("^tk_(.*?)_results.*", "\\1", colnames(tk.pval))
  
  ## Get p-values from Dunn test
  dn.tests <- list.files(dn.dir, pattern = '_processed.csv', full.names = TRUE)
  dn.pval <- as.data.frame(lapply(dn.tests, function(x) read.csv(x)$P.adj))
  names(dn.pval) <- file_path_sans_ext(basename(dn.tests))
  
  ## Fix column names for Dunn p-values
  colnames(dn.pval) <- gsub("^dn_(.*?)_results.*", "\\1", colnames(dn.pval))
  
  ## Get comparison labels from one of the processed Tukey test files
  files <- list.files(path = tk.dir, pattern = "^tk_.*_results_processed\\.csv$", full.names = TRUE)
  
  ## Ensure that at least one file is found
  if (length(files) > 0) {
    ## Read the first matching file
    df <- read.csv(files[1])
    
    ## Extract the 'Normalized_Comparison' column
    comparisons <- df$Normalized_Comparison
  } else {
    stop("No processed Tukey test files found in the specified directory.")
  }
  
  ## Combine p-values
  phoc.pval <- cbind(comparisons, tk.pval, dn.pval)
  
  ## Round all columns to 3 decimals and add asterisk if `asterisk = TRUE` and p-value < 0.05
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

  ## Only write the combined p-values to the output file if `output.file` is not NULL or not missing
  if (!is.null(output.file)) {
    write.csv(phoc.pval, output.file, row.names = FALSE)
    cat("Processed p-values written to:", output.file, "\n")
  } else {
    cat("Output file not specified. P-values not written.\n")
  }

  return(phoc.pval)
}
