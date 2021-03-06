% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motbf.R
\name{univMoTBF}
\alias{univMoTBF}
\title{Fitting MoTBFs}
\usage{
univMoTBF(
  data,
  POTENTIAL_TYPE,
  evalRange = NULL,
  nparam = NULL,
  maxParam = NULL
)
}
\arguments{
\item{data}{A \code{"numeric"} vector.}

\item{POTENTIAL_TYPE}{A \code{"character"} string specifying the potential
type, must be either \code{"MOP"} or \code{"MTE"}.}

\item{evalRange}{A \code{"numeric"} vector that specifies the domain over
which the model will be fitted. By default, it is \code{NULL} and the function is defined
over the complete data range.}

\item{nparam}{The exact number of basis functions to be used. By default, it is \code{NULL}
and the best MoTBF is fitted taking into account the Bayesian information
criterion (BIC) to score and select the functions. It evaluates the next two functions and,
if the BIC value does not improve, the function with the best BIC score so far is returned.}

\item{maxParam}{A \code{"numeric"} value which indicates the maximum number of coefficients in the function. 
By default, it is \code{NULL}; otherwise, the function which gets the best BIC score
with at most this number of parameters is returned.}
}
\value{
\code{univMoTBF()} returns an object of class \code{"motbf"}. This object is a list containing several elements, 
including its mathematical expression and other hidden elements related to the learning task. 
The processing time is one of the values returned by this function and it can be extracted by $Time. 
Although the learning process is always the same for a particular data sample, 
the processing can vary inasmuch as it depends on the CPU.
}
\description{
Function for fitting univariate mixture of truncated basis functions.
Least square optimization is used to minimize the quadratic 
error between the empirical cumulative distribution and the estimated one.
}
\examples{
## 1. EXAMPLE
## Data
X <- rnorm(5000)

## Learning
f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE"); f1
f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f2

## Plots
hist(X, prob = TRUE, main = "")
plot(f1, xlim = range(X), col = 1, add = TRUE)
plot(f2, xlim = range(X), col = 2, add = TRUE)

## Data test
Xtest <- rnorm(1000)
## Filtered data test
Xtest <- Xtest[Xtest>=min(X) & Xtest<=max(X)]

## Log-likelihood
sum(log(as.function(f1)(Xtest)))
sum(log(as.function(f2)(Xtest)))

## 2. EXAMPLE
## Data
X <- rchisq(5000, df = 5)

## Learning
f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", nparam = 11); f1
f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", maxParam = 10); f2

## Plots
hist(X, prob = TRUE, main = "")
plot(f1, xlim = range(X), col = 3, add = TRUE)
plot(f2, xlim = range(X), col = 4, add = TRUE)

## Data test
Xtest <- rchisq(1000, df = 5)
## Filtered data test
Xtest <- Xtest[Xtest>=min(X) & Xtest<=max(X)]

## Log-likelihood
sum(log(as.function(f1)(Xtest)))
sum(log(as.function(f2)(Xtest)))

}
