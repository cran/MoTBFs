% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{preprocessedData}
\alias{preprocessedData}
\title{Data cleaning}
\usage{
preprocessedData(data, strangeElements)
}
\arguments{
\item{data}{A dataset of class \code{"matrix"} or \code{"data.frame"},}

\item{strangeElements}{A \code{"character"} string which contains the elementes to remove.}
}
\description{
Delete rows of a dataset wich contains anomalous values.
}
