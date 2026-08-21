#' Fitting mixtures of polynomials
#' 
#' These functions fit mixtures of polynomials (MOPs).  
#' Least square optimization is used to  
#' minimize the quadratic error between the empirical 
#' cumulative distribution and the estimated one. 
#' 
#' @name mop.learning
#' @rdname mop.learning
#' @param X A \code{"numeric"} data vector.
#' @param nparam Number of parameters of the function.
#' @param domain A \code{"numeric"} containing the interval over which the function is defined.
#' @param maxParam A \code{"numeric"} value indicating the maximum number of coefficients in the function. By default it is \code{NULL},
#' which means that the number of parameter is not limited. The output is the function which gets 
#' the best BIC (with at most\code{maxParam} parameters if not NULL).
#' @return 
#' \code{mop.lerning()} returns a list of n elements:
#' \item{Function}{An \code{"motbf"} object of the \code{'mop'} subclass.}
#' \item{Subclass}{\code{'mop'}.}
#' \item{Domain}{The range where the function is defined to be a legal density function.}
#' \item{Iterations}{The number of iterations that the optimization problem takes to minimize
#' the errors.}
#' \item{Time}{The CPU time employed.}
#' 
#' \code{bestMOP()} returns a list including the polynomial function with the best BIC score, 
#' the number of parameters and an array with 
#' the BIC values of the evaluated functions.
#' @details 
#' \code{mop.learning()}:
#' The returned value \code{$Function} is the only visible element which contains the algebraic expression. 
#' Using \link{attributes} the name of the others elements are shown and also they can be extracted with \code{$}.
#' The \link{summary} of the function also shows all these elements.
#'
#' \code{bestMOP()}:
#' The first returned value \code{$bestPx} contains the output of the \code{mop.learning()} function
#' with the number of parameters which gets the best BIC values, taking into account the  
#' BIC score to penalize the functions. It evaluates the two next functions,
#' if the BIC score does not improve then the function with the last best BIC is returned.
#' 
#' @seealso \link{univMoTBF} A complete function for learning MOPs which includes extra options.
#' @importFrom quadprog solve.QP
#' @importFrom grDevices colorRampPalette terrain.colors
#' @importFrom graphics abline contour filled.contour layout lcm legend par persp plot points title
#' @importFrom methods is
#' @importFrom stats coef cov dunif ecdf integrate ks.test quantile rnorm runif sd uniroot
#' @importFrom utils data
#'
#' @examples
#' ## 1. EXAMPLE 
#' data <- rnorm(1000)
#' 
#' ## MOP with fix number of degrees
#' fx <- mop.learning(data, nparam=7, domain=range(data))
#' fx
#' hist(data, prob=TRUE, main="")
#' plot(fx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## Best MOP in terms of BIC
#' fMOP <- bestMOP(data, domain=range(data))
#' attributes(fMOP)
#' fMOP$bestPx
#' hist(data, prob=TRUE, main="")
#' plot(fMOP$bestPx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## 2. EXAMPLE
#' data <- rbeta(4000, shape1=1/2, shape2=1/2)
#' 
#' ## MOP with fix number of degrees 
#' fx <- mop.learning(data, nparam=6, domain=range(data))
#' fx
#' hist(data, prob=TRUE, main="")
#' plot(fx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## Best MOP in terms of BIC
#' fMOP <- bestMOP(data, domain=range(data), maxParam=6)
#' attributes(fMOP)
#' fMOP$bestPx
#' attributes(fMOP$bestPx)
#' hist(data, prob=TRUE, main="")
#' plot(fMOP$bestPx, col=2, xlim=range(data), add=TRUE)


#' @export
mop.learning <- function(X, nparam, domain)
{
  ## Time
  tm <- Sys.time()
  ## CDF quadratic function to be minimized
  x <- sort(as.numeric(X)); f <- ecdf(x); y <- f(x); 
  y  <- unique(y); x <- unique(x);  n <- length(x)
  
  ## N.records is equal 1
  if(length(x)==1){
    if(nparam!=1) return(NULL)
    P <- asMOPString(1/diff(domain))
    P <- list(Function = P, Subclass = "mop", Domain = domain,
              Iterations = 0, Time = 0)
    P <- new_mop(P)
    return(P)
  }
  xx=rep(1,n); xi <- 0
  for(i in 1:nparam){
    xi=cbind(x^i)
    xx=cbind(xx,xi)
  }
  
  ## Dmat
  XX <- t(xx)%*%xx 
  
  ## dvec
  Xy <- t(xx)%*%y
  
  ## Compute the constraints
  ma2=c()
  xnew=seq(min(domain), max(domain), 10^(-3))
  if(xnew[length(xnew)]!=max(domain)) xnew=c(xnew,max(domain))
  xi=c()
  ma=0 
  for(j in 1:nparam){
    xi=j*xnew^(j-1)
    ma=rbind(ma,xi)
  }
  ma2 <- matrix(c(ma2,ma), nrow=nparam+1)
  
  xi <- 0; ma3 <- 1
  for(i in 1:nparam){
    xi <- min(domain)^i
    ma3 <- rbind(ma3,xi)
  }
  
  xi <- 0; ma4 <- 1
  for(i in 1:nparam){
    xi <- cbind((max(domain))^i)
    ma4 <- rbind(ma4,xi)
  }
  
  ## Amat
  AA <- cbind(ma3, ma4, ma2)
  
  ## bvec
  B <- c(0, 1, rep(10^(-3),ncol(ma2)))
  
  ## Solve the optimization problem
  tr=tryCatch(solve.QP(XX, Xy, AA, B, meq=2), error = function(e) NULL)
  finaltm <- Sys.time() - tm
  if(is.null(tr)==TRUE){
    return(NULL)
  }else{
    soluc <-tr
    parameters <- soluc$solution
    Px <- asMOPString(parameters)
    Px <- list(Function=Px, Subclass="mop")
    Px <- new_mop(Px)
    
    ## Derivative (PDF)
    Px <-derivMoTBF(Px)
    
    # Garantizar que integra 1
    # browser()
    k = integrateMOP(Px,lower = domain[1],domain[2])
    
    if(abs(1-k)>10^-6){
      Px = multiplyMOPbyConstant(Px,1/k)
    }
    
    #return(Px)
    P <- asMOPString(coef(Px))
    P <- list(Function = P, Subclass = "mop", Domain = domain,
              Iterations = tr$iterations[1], Time = finaltm)
    P <- new_mop(P)
    return(P)
  }
}

 
#' @rdname mop.learning
#' @export
bestMOP=function(X, domain, maxParam=NULL)
{
  # browser()
  if(is.null(maxParam)){
    maxParam = 10
  } 
  bestBIC <- -10^10;  bestPx <- 0; degreeN <- 0
  degree <- 1; i <- 0; vecBIC <- c()
  if(!is.null(maxParam)){
    # Pxs=list()
    nparam = maxParam
    continue = TRUE
    bestPx = mop.learning(X,1,domain)
    bestBIC = sum(log(as.function(bestPx)(X)))
    while(nparam>1){
      # Control
      if(!continue){
        break
      }
      continue = FALSE
      # Learning Mop
      Px=mop.learning(X, nparam, domain)
      if(is.null(Px)) {
        nparam <- nparam-1
        continue = TRUE
        next
      }
      BiC = tryCatch(BICMoTBF(Px,X),
        # Esto puede generar NaN},
        warning = function(w) {
          continue<-TRUE
          return(-10^10)},
        error = function(e) {
          continue<-TRUE
          -10^-10
        })
      if(BiC>bestBIC){
        bestBIC = BiC
        bestPx = Px
        continue = TRUE
      }
      nparam = nparam-1
    }
    degreeN <- length(coef(bestPx))    
  } else{
    repeat{
      if(!is.null(maxParam)&&(degree>maxParam)) break
      if(i==4) break
      
      ## Learning parameters with a fix number of degree
      Px=mop.learning(X, degree, domain) #fit function
      if((is.null(Px)==TRUE)&&(!is.motbf(bestPx))) {
        degree <- degree+1
        i <- i+1
        next
      }
      if((is.null(Px)==TRUE)&&(is.motbf(bestPx))) break
      
      ## compute the BIC score
      BiC <- BICMoTBF(Px,X)
      vecBIC <- c(vecBIC,BiC)
      if(is.na(BiC)) break
      
      if(length(vecBIC)<=2){
        if(BiC>bestBIC){
          bestBIC <- BiC
          bestPx <- Px
          degreeN <- degree-1
          degree <- degree+1
        } else {
          degree <- degree+1
          next
        }
      } else{
        if(BiC>bestBIC){
          bestBIC <- BiC
          bestPx <- Px
          degreeN <- degree-1
          degree <- degree+1
        } else{
          if(bestBIC==vecBIC[length(vecBIC)-1]){
            if(length(unique(vecBIC))==1) break
            degree <- degree+1
            next
          } else{
            break
          }
        }
      }
    }
    
  }
  
  result <- list(bestPx=bestPx, bestBIC=bestBIC, degreeN=degreeN-1, vecBIC=vecBIC)
  return(result)
}


# bestMOP=function(X, domain, maxParam=NULL)
# {
#   bestBIC <- -10^10;  bestPx <- 0; degreeN <- 0
#   degree <- 1; i <- 0; vecBIC <- c()
#   if(!is.null(maxParam)){
#     Pxs=list()
#     repeat{
#       if(!is.null(maxParam)&&(degree>maxParam)) break
#       if(i==4) break
#       
#       ## Learning parameters with a fix number of degree
#       
#       Px=mop.learning(X, degree, domain) #fit function
#       if((is.null(Px)==TRUE)&&(!is.motbf(bestPx))) {
#         degree <- degree+1
#         i <- i+1
#         next
#       }
#       if((is.null(Px)==TRUE)&&(is.motbf(bestPx))) break
#       
#       ## compute the BIC score
#       degree <- degree+1
#       
#       BiC <- sum(log(as.function(Px)(X)))
#       #BiC <- BICMoTBF(Px,X)
#       Pxs[[length(Pxs)+1]] <- Px
#       vecBIC <- c(vecBIC,BiC)
#     }
#     p <- which(vecBIC==max(vecBIC))
#     bestBIC <- vecBIC[p]
#     bestPx <- Pxs[[p]]
#     degreeN <- length(coef(bestPx))    
#   } else{
#     repeat{
#       if(!is.null(maxParam)&&(degree>maxParam)) break
#       if(i==4) break
#       
#       ## Learning parameters with a fix number of degree
#       Px=mop.learning(X, degree, domain) #fit function
#       if((is.null(Px)==TRUE)&&(!is.motbf(bestPx))) {
#         degree <- degree+1
#         i <- i+1
#         next
#       }
#       if((is.null(Px)==TRUE)&&(is.motbf(bestPx))) break
#       
#       ## compute the BIC score
#       BiC <- BICMoTBF(Px,X)
#       vecBIC <- c(vecBIC,BiC)
#       if(is.na(BiC)) break
#       
#       if(length(vecBIC)<=2){
#         if(BiC>bestBIC){
#           bestBIC <- BiC
#           bestPx <- Px
#           degreeN <- degree-1
#           degree <- degree+1
#         } else {
#           degree <- degree+1
#           next
#         }
#       } else{
#         if(BiC>bestBIC){
#           bestBIC <- BiC
#           bestPx <- Px
#           degreeN <- degree-1
#           degree <- degree+1
#         } else{
#           if(bestBIC==vecBIC[length(vecBIC)-1]){
#             if(length(unique(vecBIC))==1) break
#             degree <- degree+1
#             next
#           } else{
#             break
#           }
#         }
#       }
#     }
#     
#   }
#   
#   result <- list(bestPx=bestPx, bestBIC=bestBIC, degreeN=degreeN+1, vecBIC=vecBIC)
#   return(result)
# }

#' Parameters to MOP String
#' 
#' This function builds a string with the structure of a \code{'mop'} function.
#' 
#' @param parameters A \code{"numeric"} vector containing the coefficients.
#' @return A \code{"character"} string with a \code{'mop'} structure.
#' @export
#' @examples
#' 
#' param <- c(1,2,3,4)
#' asMOPString(param)
#' 
#' param <- 3.4
#' asMOPString(param)
#' 
asMOPString <- function(parameters) 
{
  str <- parameters[1]
  if(length(parameters)==1) return(noquote(paste(str,"+0*x", sep="")))
  sign <- parameters; sign[sign<0] <- ""; sign[sign>=0] <- "+"
  for(i in 2: length(parameters)) str <- noquote(paste(str,sign[i], parameters[i], "*x",ifelse((i-1)!=1,paste("^", i-1, sep=""),""), sep=""))
  return(str)  
}


#' Extract coefficients from MOPs
#' 
#' It extracts the parameters of the learned mixtures of polynomial models.
#' 
#' @name coef.mop
#' @rdname coef.mop
#' @param fx An \code{"motbf"} function of subclass \code{'mop'}.
#' @return An array with the parameters of the function.
#' @details
#' \code{coeffMOP()} return the coefficients of the terms in the function.
#' 
#' \code{coeffPol()} returns the coefficients of the potential of the polynomial basis in the function. 
#' @seealso \link{coef.motbf} and \link{univMoTBF}
#' @examples
#' 
#' ## 1. EXAMPLE
#' data <- rchisq(1000, df=5)
#' fx1 <- univMoTBF(data, POTENTIAL_TYPE = "MOP")
#' hist(data, prob=TRUE, main="")
#' plot(fx1, xlim=range(data), col="red", add=TRUE)
#' coeffMOP(fx1) ## coef(fx1)
#' coeffPol(fx1)
#' 
#' ## 2. EXAMPLE
#' data <- rexp(1000, rate=1/2)
#' fx2 <- univMoTBF(data, POTENTIAL_TYPE = "MOP")
#' hist(data, prob=TRUE, main="")
#' plot(fx2, xlim=range(data), col="red", add=TRUE)
#' coeffMOP(fx2) ## coef(fx2)
#' coeffPol(fx2)



#' @export
coeffMOP <- function(fx)
{
  fx <- noquote(as.character(fx))
  #var <- nVariables(fx)
  mu <- tryCatch(meanMOP(fx), error = function(e) NA) 
  if(!is.na(mu)){
    fx <-strsplit(fx, split=mu)[[1]]
  }
  
  f1 <- substr(fx[1], 1, 1)
  t <- unlist(sapply(1:length(fx), function(i) strsplit(fx[i], split="-", fixed = TRUE)[[1]]))
  for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")
  if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
  
  t2 <- c()
  for(i in 1:(length(t))){
    t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  
  t3 <- c()
  for(i in 1:(length(t2))){
    t1 <- strsplit(t2[i], split="*", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
    t3 <- c(t3,t1[1])
  }
  
  pos <- grep("e", t3)
  if(length(pos)!=0){
    h1 <- c()
    for(i in pos){
      if(as.numeric(t3[i+1])>=0) sign="+" else sign=""
      h=paste(t3[i], sign, t3[i+1], sep="")
      h1=c(h1,h)
    }
    t3[pos] <- h1; t4 <- t3[-(pos+1)]; t3 <- t4
  }
  
  t3 <- t3[!sapply(1:length(t3), function(i) t3[i]=="-")]
  t3 <- t3[!sapply(1:length(t3), function(i) t3[i]=="+")]
  coeff <- as.numeric(t3)[!sapply(as.numeric(t3), is.na)]
  return(coeff)
}


#' @rdname coef.mop
#' @export
coeffPol <- function(fx)
{
  param <- coef(fx)
  string <- noquote(as.character(fx))
  for(i in 1:length(param)) string <- unlist(strsplit(string, split=param[i], fixed=TRUE))
  string <- unlist(strsplit(string, split="+", fixed=TRUE))
  if(string[1]=="") string <- string[-1]
  if(length(string)==(length(param)-1)) coeff <- 0
  else coeff <- c()
  string <- strsplit(string, split="^", fixed=TRUE)
  t <- sapply(1:length(string), function(i) string[[i]][2])
  t[is.na(t)] <- 1
  coeff <- as.numeric(c(coeff, t))
  return(coeff)
}


#' Derivative of a MOP
#' 
#' Compute the derivative of an \code{"motbf"} object with \code{'mop'} subclass.
#' 
#' @param fx An \code{"motbf"} object of the \code{'mop'} subclass.
#' @return The derivative which is also an \code{"motbf"} function.
#' @seealso \link{univMoTBF} for learning and \link{derivMoTBF} for 
#' general \code{"motbf"} models.
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMOP(Px)
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMOP(Px)
#' 
#' \dontrun{
#' ## 3. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' derivMOP(Px)
#' ## Error in derivMOP(Px): fx is an 'motbf' function but not 'mop' subclass.
#' class(Px)
#' subclass(Px)
#' }

derivMOP <- function(fx) 
{
  if(!is.motbf(fx)) stop("fx is not an 'motbf' function.")
  if(is.motbf(fx)&&!is.mop(fx)) stop("fx is an 'motbf' function but not 'mop' subclass.")
  
  parameters <- coef(fx)
  str <- parameters[2]
  if(length(parameters)<3) str <- paste(str,"+0*x", sep="")
  if(length(parameters)>=3){
    sign <- parameters; sign[sign<0]=""; sign[sign>=0]="+"
    for(i in 3: length(parameters)) str <- noquote(paste(str,sign[i], parameters[i]*(i-1), "*x",ifelse((i-2)!=1,paste("^", i-2, sep=""),""), sep=""))
  }
  str <- noquote(str)
  str <- list(Function = str, Subclass= "mop")
  str <- new_mop(str)
  return(str)  
}


#' Expected Value of an MoP Density Function
#'
#' Computes the expected value (mean) of a Mixture of Polynomials (MoP) function over its domain.
#'
#' @param fx An object of class \code{'mop'} representing an MoP probability density function.
#'
#' @return A numeric value representing the expected value of the distribution.
#' @export
#' 
expectedValueMOP = function(fx){
  f = fx
  if('piecewisemop' %in% class(f)){
    class(f) = 'piecewisemop'
    res = sapply(f, expectedValueMOP)
    return(sum(res))
  }
  # get coefficients and exponents
  des = splitMOP(f)
  coef = des[[1]]
  expo = des[[2]]
  
  # variable
  v = des[[3]]
  
  # get degree of each term
  b = ifelse(expo == "", 0,expo)# exponente para término independiente
  b = gsub(paste0('\\*',v, '|\\^'),'',b)
  b = as.numeric(ifelse(b == '', 1, b))
  
  
  lower = f$Domain[1]
  upper = f$Domain[2]
  
  fx <- function(x){sum(coef/(b+2)*x^(b+2))}
  return(fx(upper) - fx(lower))
}


expectedValue = function(fx){
  f = fx
  if('piecewisemop' %in% class(f)){
    class(f) = 'piecewisemop'
    res = sapply(f, expectedValue)
    return(sum(res))
  }
  l = f$Domain[1]
  u = f$Domain[2]
  
  f = as.function(f)
  g = function(x){ return(x*f(x))};
  
  
  res = integrate(g,l,u)
  return(res$value);
  
}

splitMOP <- function(mop){
  
  # Extraer bases y exponentes
  coefmop = coeffMOP(mop)
  str = as.character(abs(coefmop));str
  
  test = tryCatch({as.character(mop$Function)}, error=function(e){NULL})
  if(is.null(test) & is.numeric(mop) & length(mop)==1){
    test = mop
  }

  
  # remove sci notation
  test = gsub('[0-9]e-[0-9]','', test)
  test = gsub('[0-9]e\\+[0-9]','', test)
  
  # Split each term by + or -
  terms = strsplit(test, split = c("(-)|(\\+)"))[[1]]
  # If first coef is negative, remove first empty element
  if(coefmop[1]<0){
    terms = terms[-1]
  }
  # split by * within each term
  
  terms = strsplit(terms, split = '(?<=.)(?=\\*)', perl = TRUE)
  
  # keep exponents (if mop is a univ distribution, they are in the 2nd element; 
  # if mop is s joint distribution, they are in the 2nd to n element: therefore, use grep())
  # if("univmotbf"%in%class(mop)){
  #   exponents = sapply(terms, '[',2);exponents
  #   exponents = ifelse(is.na(exponents), '', exponents);exponents
  # }else if("jointmotbf"%in%class(mop)){
  #   exponents = sapply(terms, grep, pattern = "\\*", value = TRUE)
  #   exponents = lapply(exponents, function(x) ifelse(length(x)==0, '', x))
  #   exponents = sort(unique(unlist(exponents)))
  # }
  exponents = sapply(terms, grep, pattern = "\\*", value = TRUE)
  exponents = lapply(exponents, function(x) if(length(x)==0){''}else{x})
  expo = sapply(exponents, paste0, collapse = '')

  # get variables
  variables = sort(unique(unlist(exponents)))
  lp = unlist(variables)
  
  tm1 = strsplit(lp, '(\\*|\\^)')
  
  tm2 = suppressWarnings({sapply(tm1, function(x)ifelse(!is.na(as.numeric(x)),"", x))})
  
  tm3 = unique(unlist(tm2))
  res = tm3[tm3!=""]
  vars = unique(res)
  if(length(vars)==0){
    vars = ""
  }
  return(list(Coefficients = coefmop, Exponents = expo, Variables = vars )) 
}
