% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint.R
\name{as.function.jointmotbf}
\alias{as.function.jointmotbf}
\title{Coerce a \code{"jointmotbf"} Object to a Function}
\usage{
\method{as.function}{jointmotbf}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"joinmotbf"}.}

\item{\dots}{Further arguments to be passed to or from the method. Not necessary for this method.}
}
\value{
It returns a function to evaluate an object of class \code{"jointmotbf"}.
}
\description{
Takes a \code{"jointmotbf"} object and contructs an \R function to evaluate it at multidimensional points.
}
\details{
This is an \code{S3} method for the generic function \link{as.function}.
}
\examples{
## 1.EXAMPLE
## Dataset
data <- data.frame(X = rnorm(100), Y = rexp(100))

## Joint function
dim <- c(3,2)
param <- parametersJointMoTBF(data, dimensions = dim)
P <- jointMoTBF(param)
density <- as.function(P)(data[,1], data[,2])
density

## Log-likelihood
sum(log(density))

#############################################################################
## MORE EXAMPLES ############################################################
#############################################################################
\donttest{
## Dataset
data <- data.frame(X = rnorm(100), Y = rexp(100), Z = rnorm(100))

## Joint function
dim <- c(2,3,4)
param <- parametersJointMoTBF(data, dimensions = dim)
P <- jointMoTBF(param)
density <- as.function(P)(data[,1], data[,2], data[,3])
density

## Log-likelihood
sum(log(density))
}
}
\seealso{
\link{parametersJointMoTBF} and \link{jointMoTBF}
}
