% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fda.R
\name{plot_fda}
\alias{plot_fda}
\title{Plot Flexible Discriminant Analysis}
\usage{
plot_fda(
  data,
  fixed.shape = NULL,
  point.color = NULL,
  split = NULL,
  accu = TRUE
)
}
\arguments{
\item{data}{A data frame with population or species label in the first column, followed by morpholigcal data.}

\item{fixed.shape}{A vector specifying the point shape for each population or species. If NULL, random points will be used.}

\item{point.color}{A vector specifying the point color for each population or species. If NULL, random colors will be used}

\item{split}{Numeric. Proportion of data to be assigned as training data. If NULL, all data will be used as training set.}

\item{accu}{Logical. If TRUE, it prints the accuracy of classification based on confusion matrix.}
}
\value{
A biplot showing the results of FDA.
}
\description{
A function to plot the first two dimensions of Flexible Discriminant Analaysis (FDA).
}
\examples{
data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))
data$Sex <- NULL

point.shape <- c("Luzon" = 8, "Mindanao" = 11, "Palawan" = 10)
point.color <- c("Luzon" = "#000000", "Mindanao" = "#FF7F0E", "Palawan" = "#D62728")

fda.plot <- plot_fda(data, point.color = point.color, 
                               fixed.shape = point.shape, split = 0.7, accu = TRUE)
fda.plot

}
