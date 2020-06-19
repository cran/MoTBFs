#' Random generation for MoTBF distributions
#' 
#' Random generation for mixtures of truncated basis functions defined in a specific domain.
#' The inverse transform method is used.
#' 
#' @name MoTBF-Distribution
#' @rdname MoTBF-Distribution
#' @param size A non-negative integer indicating the number of records to generate.
#' @param fx An object of class \code{"motbf"}.
#' @param domain 
#' A \code{"numeric"} vector indicating the lower and upper limits to sample from.
#' If not specified, the range is taken from the object \code{fx}.
#' @param data A \code{"numeric"} vector to be compared with the simulated sample. 
#' By default, it is \code{NULL}; otherwise, the empirical cumulative distributions of both 
#' the data and the simulated sample are plotted and the Kolmogorov Smirnov test 
#' is used to test whether or not both samples can be considered to be drawn from the same distribution.
#' @return \code{rMoTBF()} returns a \code{"numeric"} vector containing the simulated values. 
#' \code{inversionMethod()} returns a list with the simulated values and the results 
#' of the two-sample Kolmogorov-Smirnov test, as well as the plot of the CDFs of the 
#' original and simulated data.
#' @seealso \link{integralMoTBF}
#' @examples
#' 
#' ## 1. EXAMPLE
#' ## Data
#' X <- rnorm(1000, mean = 5, sd = 3)
#' 
#' ## Learning
#' f <- univMoTBF(X, POTENTIAL_TYPE="MOP", nparam=10)
#' plot(f, xlim = f$Domain)
#' 
#' ## Random sample
#' Y <- rMoTBF(size = 500, fx = f)
#' ks.test(X,Y)
#' 
#' ## Plots
#' hist(Y, prob = TRUE, add = TRUE)
#' 
#' ## 2. EXAMPLE 
#' ## Data
#' X <- rweibull(5000, shape=2)
#' 
#' ## Learning
#' f <- univMoTBF(X, POTENTIAL_TYPE="MOP", nparam=10)
#' plot(f, xlim = f$Domain)
#' 
#' ## Random sample
#' inv <- inversionMethod(size = 500, fx = f, data = X)
#' attributes(inv)
#' inv$test
#' Y <- inv$sample 
#' 
#' ## Plots
#' plot(f, xlim = f$Domain)
#' hist(Y, prob = TRUE, add = TRUE)
#' 

#' @rdname MoTBF-Distribution
#' @export
rMoTBF <- function(size, fx, domain = NULL)
{
  if(is.null(domain)) domain <- fx$Domain
  if(is.null(domain)) stop("Domain is required.")
  CDF <- integralMoTBF(fx)
  intmin <- as.function(CDF)(min(domain))
  mUnif <- runif(size)
  simulatedValues <- sampleMoTBF <- sapply(1:size, function(i)
                     uniroot(as.function(motbf(paste(CDF,ifelse
                     ((-intmin-mUnif[i])>=0, "+", ""),-intmin-
                     mUnif[i], sep=""))), range(domain))$root)
  return(simulatedValues)
}

#' @name MoTBF-Distribution
#' @export
inversionMethod <- function(size, fx, domain = NULL, data = NULL)
{
  opar <- par(no.readonly =TRUE)       
  on.exit(par(opar))     
  
  simulatedValues <- rMoTBF(size, fx, domain)
  
  plot(ecdf(simulatedValues), cex = 0, main = "")
  if(!is.null(data)){
    plot(ecdf(data), col="red", cex = 0, main = "", add = T)
    par(xpd = TRUE)
    legend(0, 1.4, c(expression(F(X)),expression(F(Simulated_Values))), 
         cex=0.8, col = c("red", "black"), lty = c(1,1), lwd = c(1,1), 
         bty = "n")
    test <- ks.test(data, simulatedValues)
    simulatedValues <- list(sample = simulatedValues, test = test)
  }
  
  return(simulatedValues)
}

