% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_sex.R
\name{filter_sex}
\alias{filter_sex}
\title{Filter Data by Sex}
\usage{
filter_sex(data)
}
\arguments{
\item{data}{A dataframe with population and sex information in the first two columns, followed by morphology data.}
}
\value{
A data frame containing filtered data.
}
\description{
A function to filter data by excluding populations represented by only one sex.
Additionally, it filters sex data within each population with fewer than three samples.
}
\examples{
data <- read.csv(system.file("extdata", "herp.data.csv", package = "MorphoTax"))

filtered.data <- filter_sex(data) 

filtered.data

}
