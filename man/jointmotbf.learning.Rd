% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/joint.R
\name{jointmotbf.learning}
\alias{jointMoTBF}
\alias{jointmotbf.learning}
\alias{parametersJointMoTBF}
\title{Learning Joint Functions}
\usage{
parametersJointMoTBF(X, ranges = NULL, dimensions = NULL)

jointMoTBF(object)
}
\arguments{
\item{X}{A dataset of class \code{"data.frame"}.}

\item{ranges}{A \code{"numeric"} matrix by columns with the ranges of the varibles where fitting the function.}

\item{dimensions}{A \code{"numeric"} vector with the number of parameters of each varible.}

\item{object}{A list with the output of the function \code{parametersJointMoTBF()}.}
}
\value{
\code{parametersJointMoTBF()} returns a list with the elements: \bold{Parameters}, which contains the computed
coefficients of the resulting function;
\bold{Dimension} which is a \code{"numeric"} vector containing the number
of coefficients used for each variable;
\bold{Range} contains a \code{"numeric"} matrix with the domain of each variable by columns;
\bold{Iterations} contains a number of iterations needed to solve the problem;
\bold{Time} contains the time that the functions spent to solve the problem.

\code{jointMoTBF()} retunrs an object of class \code{"jointmotbf"}; It is a list whose unique visible element
is the mathematical expression, furthermore it contains the other elements of the output of the
\code{parametersJointMoTBF()} function.
}
\description{
Two functions for learning joint MoTBFs. The first one, \code{parametersJointMoTBF()},
gets the parameters by solving a quadratic optimization problem, minimizing
the error of the empirical cumulative joint values versus the estimated ones.
Finally it derives these parameters in order to obtain parameters of the joint density function.
The second one, \code{jointMoTBF()}, fixes the equation of the joint function using
the previously learned parameters and converting this \code{"character"} string in an
object of class \code{"jointmotbf"}.
}
\examples{
## 1. EXAMPLE
require(mnormt) ## Load the package to generate a multinormal dataset
## Dataset
Mean <- 0; nVar <- 2; ro <- 0
varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
means <- rep(Mean, nVar)
X <- rmnorm(100,means,varcov)
data <- standardizeDataset(data.frame(X))

## Joint learnings
dim <- c(2,3)
param <- parametersJointMoTBF(X = data, dimensions = dim)

param$Parameters
length(param$Parameters) ## prod(dim)
param$Dimension
param$Range

P <- jointMoTBF(param)
P
attributes(P)
class(P)

###############################################################################
## MORE EXAMPLES ##############################################################
###############################################################################
## ## Dataset
## Mean <- 1; nVar <- 3; ro <- 0.5
## varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
## means <- rep(Mean, nVar)
## X <- rmnorm(200,means,varcov)
## data <- standardizeDataset(data.frame(X))
##
## ## Joint learnings
## dim <- c(3,2,4,2)
## param <- parametersJointMoTBF(X = data, dimensions = dim)
##
## param$Parameters
## length(param$Parameters) ## prod(dim)
## param$Dimension
## param$Range
## param$Time
##
## P <- jointMoTBF(param)
## P
## attributes(P)
## class(P)
###############################################################################
###############################################################################
}
