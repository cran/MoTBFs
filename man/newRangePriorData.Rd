% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/priorKnowledge.R
\name{newRangePriorData}
\alias{newRangePriorData}
\title{Redefining the Domain}
\usage{
newRangePriorData(fPI, priorData, N, domain, s, POTENTIAL_TYPE)
}
\arguments{
\item{fPI}{The fitted function to the prior data of class \code{"motbf"}.}

\item{priorData}{A \code{"numeric"} array with the values we want to include as prior information.}

\item{N}{A \code{"numeric"} value which is the size of the data.}

\item{domain}{A \code{"numeric"} array with the limits where defining the data function.}

\item{s}{A \code{"numeric"} value which is the confident of the expert in his information. It is between 0 and the data size.}

\item{POTENTIAL_TYPE}{A \code{"character"} string giving the potential of the model, i.e. \code{"MOP"} if the basis functions are polynomials,
or \code{"MTE"} if they are exponentials.}
}
\value{
A \code{"numeric"} array which contains the new domain of the prior function.
}
\description{
Computes the new domain of two datasets.
}
\examples{
## Data
X <- rnorm(15)

## Prior Data
priordata <- rnorm(5000)

## Learning
type = "MTE" ## type <- "MOP"
fPrior <- univMoTBF(priordata, POTENTIAL_TYPE = type)

## New range
confident <- 5 ## confident <- 1,2,...,length(X)
domain <- range(X)
N <- length(X)
newRange <- newRangePriorData(fPrior, priorData = priordata, N = N,
domain = domain, s = confident, POTENTIAL_TYPE = type)
newRange
}
