#' Integrating MoTBFs
#' 
#' Compute the integral of a one-dimensional mixture of truncated basis function 
#' (objects of class \code{"mop"} or \code{"mte"}) over a bounded or unbounded interval, 
#' or compute the indefinite integral of a joint function  (object of class \code{"jointmotbf"})
#' over a subset of variables or over all the variables in the function.
#' 
#' @param f An object of class \code{"mop"}, \code{"mte"} or \code{"jointmotbf"}.
#' @param ... optional arguments to be passed to subsequent methods for the integrate.motfb() function. See details.
#' @details 
#' This function is a wrapper of internal functions integrateMOP(), integrateMTE() and integrateJointmotbf().
#' 
#' If f is of class \code{"mop"} or \code{"mte"}, valid optional arguments are 'lower' 
#' (the lower integration limit) and 'upper' (the upper integration limit), 
#' which represent the limits of the interval to compute the definite integral. 
#' If 'lower' and 'upper' are not specified, then the output is the expression of the indefinite integral.
#' 
#' On the other hand, if f is of class \code{"jointmotbf"}, the only valid optional argument is 'var', which
#' is a \code{"character"} vector containing the name of the variables that will be integrated out.
#' If not specified, then all the variables are integrated out.
#' 
#' Note that integrate.motfb() deprecates the following functions, included in previous versions 
#' of the package: \link{integralMTE}, \link{integralMOP}, \link{integralMoTBF} and \link{integralJointMoTBF}.
#' 
#' @return If f is of class \code{"mop"} or \code{"mte"}, \code{integrate.motfb()} returns 
#' either the indefinite integral of the MoTBF function, which is also an object of classes \code{"motbf"},  
#' \code{"univmotbf"}, and either \code{"mop"} or  \code{"mte"}); or the definite integral, which is a \code{"numeric"} value.
#' If f is of class \code{"jointmotbf"}, \code{integrate.motfb()} returns a multi-integral of the joint function,
#' which is also of class \code{"jointmotbf"}.
#' @seealso \link{univMoTBF} and \link{jointMoTBF}
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' ## Univariate MOP integral
#' X <- rexp(1000)
#' fx <- univMoTBF(X, POTENTIAL_TYPE = "MOP", scale = FALSE)
#' integrate.motbf(fx)
#' integrate.motbf(fx, 2, 6)
#' 
#' 
#' ## 2. EXAMPLE
#' ## Univariate MOP integral and plot of result
#' Y <- rnorm(1000)
#' fy <- univMoTBF(Y, POTENTIAL_TYPE = "MOP")
#' Fy <- integrate.motbf(fy)
#' plot(Fy)
#' integrate.motbf(fy, min(Y), max(Y))
#' 
#' ## 3. EXAMPLE
#' ## Univariate MTE integral
#' Z <- rchisq(1000, df = 3)
#' fz <- univMoTBF(Z, POTENTIAL_TYPE = "MTE", scale = FALSE)
#' integrate.motbf(fz)
#' integrate.motbf(fz, lower = 2, upper = 5)
#' 
#' \dontrun{
#' ## 4. EXAMPLE
#' Px <- "1+x+5"
#' class(Px)
#' integrate.motbf(Px)
#' ## Error in integrate.motbf(Px) : Argument "f" is not of class "motbf"
#'}
#'
#'## 5. EXAMPLE: Joint MOP integral
#'## Dataset with 2 variables
#'data <- data.frame(x = rnorm(100), y = rnorm(100))
#'
#'## Joint function
#'dim <- c(2, 3)
#'P <- jointmotbf.fit(data, dimensions = dim)
#'
#'## Integral
#'integrate.motbf(P)
#'integrate.motbf(P, var = "x")
#'integrate.motbf(P, var = "y")
#'
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' \donttest{
#' ## Dataset with 3 variables
#' data <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
#' 
#' ## Joint function
#' dim <- c(2,2,3)
#' P = jointmotbf.fit(data, dimensions = dim)
#'  
#' ## Integral
#' integrate.motbf(P)
#' integrate.motbf(P, var="x")
#' integrate.motbf(P, var=c("x","z"))
#' }
#' 
integrate.motbf <- function(f, ...){
  if(is.jointmotbf(f)){
    integrateJointmotbf(f, ...)
  }else if(is.mop(f)){
    integrateMOP(f, ...)
  }else if(is.mte(f)){
    integrateMTE(f, ...)
  }else{
    stop('Argument "f" is not of class "motbf"')
  }
}





# Integrating MOPs -------
#' Integration of MOPs 
#' 
#' Method to calculate the defined or non-defined integral of an \code{"motbf"} object of \code{'mop'} subclass.
#' 
#' @param f An \code{"motbf"} object of subclass \code{'mop'}.
#' @param lower the lower integration limit for definite integrals. By default, it is NULL.
#' @param upper the upper integration limit for definite integrals. By default, it is NULL.
#' @return The defined or non-defined integral of the function.
#' @seealso \link{univMoTBF} for learning and \link{integrate.motbf} 
#' for a more complete function to get defined and non-defined integrals
#' of class \code{"motbf"}.
#' @noRd

integrateMOP = function(f, lower = NULL, upper = NULL){
  if(!is.mop(f)) stop("Argument 'f' is  not of class 'mop'")
  
  # get coefficients and exponents
  param = splitMOP(f)
  # variable
  x = param$Variables
  
  a = param$Coefficients
  b = param$Exponents
  b = ifelse(b == "", 0, b)# exponente para término independiente
  b = gsub(paste0('\\*',x, '|\\^'),'',b)
  b = as.numeric(ifelse(b == '', 1, b))
  
  if(is.null(lower)|is.null(upper)){
    
    sn = ifelse(a <0, '', '+')
    sn[1] = ''
    
    exponentes = ifelse(b+1>1, paste0(x, "^", b+1), paste0(x))
    
    # Indefinite integral
    fx = paste0(sn, a/(b+1),"*",exponentes, collapse = '')
    
    f$Function = noquote(fx)
    return(f)
  }else{
    # Definite integral
    fx <- function(x){sum(a/(b+1)*x^(b+1))}
    return(fx(upper) - fx(lower))
  }
}



# Integrating MTEs -------
#' Integrating MTEs
#' 
#' Method to calculate the defined or non-defined integral of an \code{"motbf"} object of \code{'mte'} 
#' subclass.
#' 
#' @param f An \code{"motbf"} object of subclass \code{'mte'}.
#' @param lower the lower integration limit for definite integrals. By default, it is NULL.
#' @param upper the upper integration limit for definite integrals. By default, it is NULL.
#' @return The defined or non-defined integral of the function.
#' @seealso \link{univMoTBF} for learning and \link{integrate.motbf} 
#' for a more complete function to get defined and non-defined integrals
#' of class \code{"motbf"}.
#' @noRd
integrateMTE <- function(f, lower = NULL, upper = NULL){
  if(!is.mte(f)) stop("Argument 'f' is  not of class 'mte'")
  
  param = splitMTE(f)
  x = param$Variables
  
  a0 = param$Coefficients[1]
  a = param$Coefficients[-1]
  
  b = param$Exponents
  b = gsub(paste0('\\*',x, '|\\^'),'',b);b
  b = as.numeric(ifelse(b == '', 1, b))[-1];b
  
  if(is.null(lower)|is.null(upper)){
    
    sn = ifelse(a/b <0, '', '+')
    
    # Indefinite integral
    fx = paste0(c(a0,"*", x, paste0(sn, a/b,"*exp(", b,"*",x,")")), collapse = ''); fx
    
    f$Function = noquote(fx)
    return(f)
  }else{
    
    # Definite integral
    fx <- function(x){a0*x+ sum(a/b*exp(b*x))}
    
    return(fx(upper) - fx(lower))
  }
}

# integralMTE <-  integrateMTE

# Integrating Joint distributions -------
#' Integration with MoTBFs
#' 
#' Integrate a \code{"jointmotbf"} object over an non defined domain. It is able to
#' get the integral of a joint function over a set of variables or over all
#' the variables in the function.
#' 
#' @param P A \code{"jointmotbf"} object.
#' @param var A \code{"character"} vector containing the name of the variables that will be integrated out.
#' Instead of the names, the position of the variables can be given.
#' By default it's \code{NULL} then all the variables are integrated out.
#' @return A multiintegral of a joint function of class \code{"jointmotbf"}.
#' @noRd

integrateJointmotbf <- function(P, var=NULL){
  
  if(is.null(var)){
    # variables sobre las que se integra. si no se especifica, son todas
    var <- getMotbfVar(P)
  } 
  
  # split mop in coeffient and literal part
  splitmop = splitMOP(P) 
  pl <- splitmop$Exponents # parte literal
  coefs <- splitmop$Coefficients # coefficients
  # pl1 = pl
  # coefs1 = coefs
  
  for(h in 1:length(var)){
    pl_list = strsplit(pl, '\\*') # split literal part in each term by variables
    
    for(i in 1:length(pl_list)){# for each term, 
      pos = grep(var[h], pl_list[[i]]) # find the position of variable to integrate
      
      if(length(pos)==0){ # la variable sobre la que se integra no esta en el termino i
        int_lp = paste0(paste0(pl_list[[i]],collapse = '*'), '*', var[h]) # integrated literal part
        int_coef = coefs[i]
      }else{
        p = pl_list[[i]][pos]
        v = strsplit(p, '\\^')[[1]][1]
        g = strsplit(p, '\\^')[[1]][2]
        if(is.na(as.numeric(g))){
          new.grade = 2
        }else{
          new.grade = as.numeric(g)+1
        }
        if(new.grade == 1){
          intVar = v
        }else{
          intVar = paste0(v,'^',new.grade)
        }
        int_lp = paste0(paste0(pl_list[[i]][-pos],collapse = '*'), '*', intVar)
        int_coef = coefs[i]/new.grade
      }
      pl[i] = int_lp
      coefs[i] = int_coef
    }
  }
  
  
  s = paste0(coefs, pl);s
  # paste0(coefs1, pl1)
  # A = paste0(ifelse(sign(coefs)<0, '', '+'), s, collapse = '')
  A = paste0(s[1],paste0(ifelse(sign(coefs[-1])<0, '', '+'), s[-1], collapse = ''))
  mop_str <- noquote(A)
  if(length(splitmop$Variables)>1){ 
    P <- new_jointmotbf(mop_str)
  }else{
    P <- new_mop(mop_str)
  }
  return(P)
}