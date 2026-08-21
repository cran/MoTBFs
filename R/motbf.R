#' Fitting MoTBFs
#' 
#' Function for fitting univariate mixture of truncated basis functions.
#' Least square optimization is used to minimize the quadratic 
#' error between the empirical cumulative distribution and the estimated one. 
#' 
#' @param x A \code{"numeric"} vector.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the potential
#' type, must be either \code{"MOP"} or \code{"MTE"}.
#' @param evalRange A \code{"numeric"} vector that specifies the domain over
#' which the model will be fitted. By default, it is \code{NULL} and the function is defined
#' over the complete data range.
#' @param nparam The exact number of basis functions to be used. By default, it is \code{NULL}
#' and the best MoTBF is fitted taking into account the Bayesian information
#' criterion (BIC) to score and select the functions. It evaluates the next two functions and,
#' if the BIC value does not improve, the function with the best BIC score so far is returned.
#' @param maxParam A \code{"numeric"} value which indicates the maximum number of coefficients in the function. 
#' By default, it is \code{NULL}; otherwise, the function which gets the best BIC score
#' with at most this number of parameters is returned.
#' @param scale A \code{"logical"} value indicating whether to standardize the numeric vector (x) to have mean 0 and standard deviation 1.
#' @return \code{univMoTBF()} returns an object of class \code{"motbf"}. This object is a list containing several elements, 
#' including its mathematical expression and other hidden elements related to the learning task. 
#' The processing time is one of the values returned by this function and it can be extracted by $Time. 
#' Although the learning process is always the same for a particular data sample, 
#' the processing can vary inasmuch as it depends on the CPU.
#' @export
#' @examples
#' ## 1. EXAMPLE
#' ## Data
#' X <- rnorm(5000)
#' 
#' ## Learning
#' f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE"); f1
#' f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f2
#' 
#' ## Plots
#' hist(X, prob = TRUE, main = "")
#' plot(f1, xlim = range(X), col = 1, add = TRUE)
#' plot(f2, xlim = range(X), col = 2, add = TRUE)
#' 
#' ## Data test
#' Xtest <- rnorm(1000)
#' ## Filtered data test
#' Xtest <- Xtest[Xtest>=min(X) & Xtest<=max(X)]
#' 
#' ## Log-likelihood
#' sum(log(as.function(f1)(Xtest)))
#' sum(log(as.function(f2)(Xtest)))
#' 
#' ## 2. EXAMPLE
#' ## Data
#' X <- rchisq(5000, df = 5)
#' 
#' ## Learning
#' f1 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", nparam = 11); f1
#' f2 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", maxParam = 10); f2
#' 
#' ## Plots
#' hist(X, prob = TRUE, main = "")
#' plot(f1, xlim = range(X), col = 3, add = TRUE)
#' plot(f2, xlim = range(X), col = 4, add = TRUE)
#' 
#' ## Data test
#' Xtest <- rchisq(1000, df = 5)
#' ## Filtered data test
#' Xtest <- Xtest[Xtest>=min(X) & Xtest<=max(X)]
#' 
#' ## Log-likelihood
#' sum(log(as.function(f1)(Xtest)))
#' sum(log(as.function(f2)(Xtest)))
#' 

univMoTBF <- function(x, POTENTIAL_TYPE, evalRange=NULL, nparam=NULL,  maxParam=NULL, scale = TRUE)
{
  # browser()
  m = mean(x)
  s = sd(x)
  # if length(x) == 1, s = NA
  if(is.na(s)){s = 0.1}
  
  if(scale){
    x = (x-m)/s
  }
    
  if(is.null(evalRange)){ 
    evalRange <- range(x) 
  }else{ 
    evalRange <- evalRange
    if(scale){
      evalRange = (evalRange-m)/s
    }
  }
  
  if(is.null(nparam)){
    if(POTENTIAL_TYPE=="MOP"){
      P=bestMOP(x, evalRange, maxParam=maxParam)$bestPx
    } else if(POTENTIAL_TYPE=="MTE"){
      P=bestMTE(x, evalRange, maxParam=maxParam)$bestPx
    } else{
      stop("Unknown method, please use MOP or MTE \n")
    } 
  } else {
    if(POTENTIAL_TYPE=="MOP"){
      P <- mop.learning(x, nparam, evalRange)
    } else if(POTENTIAL_TYPE=="MTE"){
      P <- mte.learning(x, nparam, evalRange)
    } else{
      stop("Unknown method, please use MOP or MTE")
    } 
  }
  
  attr(P, 'mean') = m
  attr(P, 'sd') = s
  if(POTENTIAL_TYPE=='MOP' & scale==TRUE){
    P<-rescaledMOP(P)
    
    # Garantizar que integra a 1
    k = integrate.motbf(P, P$Domain[1], P$Domain[2])
    P = multiplyMOPbyConstant(P,1/k)
  }
  if(POTENTIAL_TYPE=="MTE" & scale==TRUE){
    x = x*s+m
    P<-rescaledMTE(P,x)
  }
  # class(P) = c(class(P), tolower(POTENTIAL_TYPE))
  # add class
  
  P = do.call(paste0("new_",subclass(P)), list(P))
  

  return(P)
}



  
  
#'Computing the BIC score of an MoTBF function
#'
#'Computes the Bayesian information criterion value (BIC) of a 
#'mixture of truncated basis functions. The BIC score is the log likelihood 
#'penalized by the number of parameters of the function and the number of
#'records of the evaluated data.
#'
#'@param Px A function of class \code{"motbf"}.
#'@param X A \code{"numeric"} vector with the data to evaluate.
#'@return A \code{"numeric"} value corresponding to the BIC score.
#'@seealso \link{univMoTBF}
#'@export
#'@examples
#'
#'## Data
#'X <- rexp(10000)
#'
#'## Data test
#'Xtest <- rexp(1000)
#'Xtest <- Xtest[Xtest>=min(X) & Xtest<=max(X)]
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", nparam = 10); f1
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 11); f2
#'
#'## BIC values
#'BICMoTBF(Px = f1, X = Xtest)
#'BICMoTBF(Px = f2, X = Xtest)
#'
#'
BICMoTBF <- function(Px, X){
  # browser()
  pPx  <-  as.function(Px)(X)
  size  <-  length(coef(Px))
  BiC  <- sum(log(pPx))-(1/2*(size+1)*log(length(X)))
  return(BiC)
}


#' Extract the coefficients of an MoTBF
#' 
#' Extracts the parameters of the learned mixtures of truncated basis
#' functions. 
#' 
#' @param object An object of class \code{motbf}.
#' @param \dots other arguments.
#' @return A numeric vector with the parameters of the function.
#' @seealso \link{univMoTBF}, \link{coeffMOP} and \link{coeffMTE}
#' @export
#' @examples
#'
#'## Data
#'X <- rchisq(2000, df = 5)
#'
#'## Learning
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); f1
#'## Coefficients
#'coef(f1)
#'
#'## Learning
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 10); f2
#'## Coefficients
#'coef(f2)
#'
#'## Learning
#'f3 <- univMoTBF(X, POTENTIAL_TYPE = "MOP", nparam=10); f3
#'## Coefficients
#'coef(f3)
#'
#'## Plots
#'plot(NULL, xlim = range(X), ylim = c(0,0.2), xlab="X", ylab="density")
#'plot(f1, xlim = range(X), col = 1, add = TRUE)
#'plot(f2, xlim = range(X), col = 2, add = TRUE)
#'plot(f3, xlim = range(X), col = 3, add = TRUE)
#'hist(X, prob = TRUE, add= TRUE)
#'
coef.motbf <- function(object, ...)
{
  if(is.mop(object)) parameters <- coeffMOP(object)
  if(is.mte(object)) parameters <- coeffMTE(object)
  return(parameters)
}




#' Derivating MoTBFs
#' 
#' Compute the derivative of a one-dimensional mixture of truncated basis function.
#' 
#' @param fx An object of class \code{"motbf"}.
#' @return The derivative of the MoTBF function, which is also 
#' an object of class \code{"motbf"}.
#' @seealso \link{univMoTBF}, \link{derivMOP} and \link{derivMTE}
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMoTBF(Px)
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMoTBF(Px)
#' 
#' ## 3. EXAMPLE
#' X <- rchisq(1000, df = 3)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' derivMoTBF(Px)
#' 
#' \dontrun{
#' ## 4. EXAMPLE
#' Px <- "x+2"
#' class(Px)
#' derivMoTBF(Px)
#' ## Error in derivMoTBF(Px): "fx is not an 'motbf' function."
#' }

derivMoTBF <- function(fx)
{
  if(!is.motbf(fx)) stop("fx is not an 'motbf' function.")
  if(is.mop(fx)) f <- derivMOP(fx)
  if(is.mte(fx)) f <- derivMTE(fx)
  return(f)
}

