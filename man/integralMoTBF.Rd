% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motbf.R
\name{integralMoTBF}
\alias{integralMoTBF}
\title{Integrating MoTBFs}
\usage{
integralMoTBF(fx, min = NULL, max = NULL)
}
\arguments{
\item{fx}{An object of class \code{"motbf"}.}

\item{min}{The lower integration limit. By default it is NULL.}

\item{max}{The upper integration limit. By default it is NULL.}
}
\value{
\code{integralMoTBF()} returns either the indefinite integral of the MoTBF 
function, which is also an object of class \code{"motbf"}, or the definite integral, 
wich is a \code{"numeric"} value.
}
\description{
Compute the integral of a one-dimensional mixture of truncated basis function 
over a bounded or unbounded interval.
}
\details{
If the limits of the interval, min and max are NULL, then the output is
the expression of the indefinite integral. If only 'min' contains a numeric value,
then the expression of the integral is evaluated at this point.
}
\examples{

## 1. EXAMPLE
X <- rexp(1000)
Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
integralMoTBF(Px)
integralMoTBF(Px, 1.2)
integralMoTBF(Px, min(X), max(X))

## 2. EXAMPLE
X <- rnorm(1000)
Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
iP <- integralMoTBF(Px); iP
plot(iP, xlim=range(X))
integralMoTBF(Px, 0.2)
integralMoTBF(Px, min(X), max(X))

## 3. EXAMPLE
X <- rchisq(1000, df = 3)
Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
integralMoTBF(Px)
integralMoTBF(Px, 1)
integralMoTBF(Px, min(X), max(X))

\dontrun{
## 4. EXAMPLE
Px <- "1+x+5"
class(Px)
integralMoTBF(Px)
## Error in integralMoTBF(Px): "fx is not an 'motbf' function."
}
}
\seealso{
\link{univMoTBF}, \link{integralMOP} and \link{integralMTE}
}
