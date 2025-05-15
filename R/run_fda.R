#' @title Flexible Discriminant Analysis
#'
#' @description A function to perform Flexible Discriminant Analysis (FDA), a classification method that integrates non-linear regression with
#' linear discriminant analysis. The function uses the 'fda' function from the 'mda' package. It then writes the confusion matrix,  
#' extracts the coefficients, and calculates the accuracy based on confusion matrix.
#'
#' @param data A data frame with population or species label in the first column, followed by morpholigcal data.
#' @param con.mat Logical (TRUE or FALSE). If TRUE, it writes the confusion matrix to a csv file.
#' @param accu Logical (TRUE or FALSE). If TRUE, it calculates and prints the classification accuracy.
#' @param coef Logical (TRUE or FALSE). If TRUE, it writes the coefficient values to a csv file.
#' @param split Numeric. Proportion of data to be assigned as training data. If NULL, all data will be used as training set.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
#' data$Sex <- NULL
#'
#' fda.res <- run_fda(data, con.mat = TRUE, coef = TRUE, accu = TRUE, split = 0.7)
#' 
#' fda.res
#'
#' @importFrom utils write.table
#' @importFrom tibble rownames_to_column
#' @return A list containing the FDA results.
#' @export

run_fda <- function(data, con.mat = TRUE, accu = TRUE, coef = TRUE, split = NULL) {
  
  ## Check input data structure
  if (!"data.frame" %in% class(data)) {
    stop("Input data must be a data frame")
  }
  
  ## Set the seed for reproducibility
  set.seed(123)
  
  ## Check if split is NULL, and do not split the data if true
  if (is.null(split)) {
    train_data <- data
    test_data <- data
  } else {
    shuffled_index <- sample(1:nrow(data)) 
    train_size <- floor(split * nrow(data))
    train_index <- shuffled_index[1:train_size]
    test_index <- shuffled_index[(train_size + 1):nrow(data)] 
    train_data <- data[train_index, ]
    test_data <- data[test_index, ]
  }
  
  ## Perform FDA on the training data
  fda.model <- fda(Pop ~ ., data = train_data)
  
  ## Initialize output list
  output <- list(fda.model = fda.model)
  
  ## Generate confusion matrix and write to CSV if con.mat = TRUE
  if (con.mat) {
    ## Apply the model to the test data
    test_pred <- predict(fda.model, newdata = test_data)
    con.mat <- table(test_data$Pop, test_pred) 
    
    ## Convert confusion matrix to a matrix
    con.mat_filename <- paste0("fda_conf_mat_", deparse(substitute(data)), ".csv")
    write.table(con.mat, con.mat_filename, sep = ",", col.names = NA)
    print(paste("Confusion matrix file written:", con.mat_filename))
    output$con.mat <- con.mat
  }
  
  ## Get FDA coefficients and write to CSV if coef = TRUE
  if (coef) {
    coef_data <- as.data.frame(fda.model$fit$coefficients)
    coef_data <- rownames_to_column(coef_data, var = "var")
    coef_filename <- paste0("fda_coef_", deparse(substitute(data)), ".csv")
    write.csv(coef_data, coef_filename, row.names = FALSE)
    print(paste("Coefficient file written:", coef_filename))
    output$coef <- coef_data 
  }
  
  ## Check accuracy of classification if accu = TRUE
  accuracy <- NA
  if (accu) {
    if (is.null(split)) {
      ## Apply the model to the training data
      train_pred <- predict(fda.model, newdata = train_data)
      tbl <- table(train_data$Pop, train_pred) 
      accuracy <- sum(diag(tbl)) / sum(tbl) 
      accuracy_percent <- accuracy * 100
      print(paste("FDA Classification Accuracy on Training Data:", round(accuracy_percent, 2), "%"))
      output$accuracy <- accuracy_percent 
      output$accuracy_data <- "Training Data"
    } else {
      ## Apply the model to the test data if split is not NULL
      test_pred <- predict(fda.model, newdata = test_data)
      tbl <- table(test_data$Pop, test_pred) 
      accuracy <- sum(diag(tbl)) / sum(tbl) 
      accuracy_percent <- accuracy * 100
      print(paste("FDA Classification Accuracy on Test Data:", round(accuracy_percent, 2), "%"))
      output$accuracy <- accuracy_percent  
      output$accuracy_data <- "Test Data" 
    }
  }

  return(output)
}

