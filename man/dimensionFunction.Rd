% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint.R
\name{dimensionFunction}
\alias{dimensionFunction}
\title{Dimension of MoTBFs}
\usage{
dimensionFunction(P)
}
\arguments{
\item{P}{An object of class \code{"motbf"} and subclass 'mop' or \code{"jointmotbf"}.}
}
\value{
Dimension of the function.
}
\description{
Get the dimension of \code{"motbf"} and \code{"jointmotbf"} densities.
}
\examples{
## 1. EXAMPLE 
## Data
X <- rnorm(2000)

## Univariate function
subclass <- "MOP"
f <- univMoTBF(X, POTENTIAL_TYPE = subclass)
dimensionFunction(f)

## 2. EXAMPLE 
## Dataset with 2 variables
X <- data.frame(rnorm(100), rnorm(100))

## Joint function
dim <- c(2,3)
param <- parametersJointMoTBF(X, dimensions = dim)
P <- jointMoTBF(param)

## Dimension of the joint function
dimensionFunction(P)

}
\seealso{
\link{univMoTBF} and \link{jointMoTBF}
}
