#' Rescaling MoTBF functions
#' 
#' A collection of function to reescale an MoTBF function
#' to the original offset and scale. This is useful when data was
#' standardized previously to learning.
#' 
#' @name rescaledFunctions
#' @rdname rescaledFunctions
#' @param fx A function of class \code{"motbf"} learned from a scaled data.
#' @param data A \code{"numeric"} vector containing the original data (non standardizded).
#' @param parameters A \code{"numeric"} vector with the coefficients to create the rescaled MoTBF.
#' @param num A \code{"numeric"} value which contains the denominator of the coefficient
#' in the exponential. By default it is 5.
#' @seealso \link{univMoTBF}
#' @return An \code{"motbf"} function of the original data.
#' @examples
#' ## 1. EXAMPLE
#' X <- rchisq(1000, df = 8) ## data
#' modX <- scale(X) ## scale data
#' 
#' ## Learning
#' f <- univMoTBF(modX, POTENTIAL_TYPE = "MOP", nparam=10) 
#' plot(f, xlim = range(modX), col=2)
#' hist(modX, prob = TRUE, add = TRUE)
#' 
#' ## Rescale
#' origF <- rescaledMoTBFs(f, X) 
#' plot(origF, xlim = range(X), col=2)
#' hist(X, prob = TRUE, add = TRUE)
#' expectedValueMOP(origF) 
#' mean(X)
#' 
#' ## 2. EXAMPLE 
#' X <- rweibull(1000, shape = 20, scale= 10) ## data
#' modX <- as.numeric(scale(X)) ## scale data
#' 
#' ## Learning
#' f <- univMoTBF(modX, POTENTIAL_TYPE = "MTE", nparam = 9) 
#' plot(f, xlim = range(modX), col=2, main="")
#' hist(modX, prob = TRUE, add = TRUE)
#' 
#' ## Rescale
#' origF <- rescaledMoTBFs(f, X) 
#' plot(origF, xlim = range(X), col=2)
#' hist(X, prob = TRUE, add = TRUE)
#' expectedValueMTE(origF)
#' mean(X)


#' @export
rescaledMoTBFs <- function(fx, data = NULL)
{
  if(is.mop(fx)) f <- rescaledMOP(fx, data)
  if(is.mte(fx)) f <- rescaledMTE(fx, data)
  return(f)
}

#' @rdname rescaledFunctions
#' @export

rescaledMOP=function(fx, data = NULL){
  f = fx
  if(is.null(data)){
    #Acceder a media y desviación típic
    m <- attr(f,'mean')
    s <- attr(f,'sd')
  }else{
    m <- mean(data)
    s <- sd(data)
  }

  signm = sign(m)

  x = getMotbfVar(f)
  polbase <- paste0(-signm*abs(m),'+',1,'*',x)
  g = f
  g$Function <- polbase
  g$Domain = data.frame(g$Domain)
  
  # coefs <- splitMOP(f)$Coefficients
  coefs <- coeffMOP(f)
  grado <- length(coefs)-1
  
  coef1 <- g
  coef1$Function = coefs[1]/s
  
  pols=list()
  
  # for (i in 1:(grado)) {
  #   pol=g
  #   if(i>=2){
  #     for (j in 1:(i-1)) {
  #       pol=multiply2Polynomials(pol,g)
  #     }
  #   }
  #   pol=multiplyMOPbyConstant(pol,coefs[i+1]/s^(i+1))
  #   pols[[i]]=pol
  #   
  # }

  pol=g
  for (i in 1:(grado)) {
    if(i>=2){
      pol=multiply2Polynomials(pol,g)
    }
    pols[[i]]=pol
  }

  pols = lapply(1:length(pols), function(i){
    multiplyMOPbyConstant(pols[[i]],coefs[i+1]/s^(i+1))
  })


  pol=pols[[1]]
  
  if(length(pols)>1){
    for (i in 1:(length(pols)-1)) {
      
      pol=simplify2MOPs(pol,pols[[i+1]])
    }
  }
  
  pol=simplify2MOPs(pol,coef1)
  
  a=f$Domain[1]*s+m
  b=f$Domain[2]*s+m
  # pol$Domain=data.frame(domain=c(a,b))
  
  pol$Domain = f$Domain
  pol$Domain[1] = a
  pol$Domain[2] = b
  # pol$Domain=c(a,b)
  return(pol)
  
}



#' @rdname rescaledFunctions
#' @export
rescaledMTE <- function(fx, data)
{
  # browser()
  parameters <- coef(fx)
  parExp <- coeffExp(fx)[-1]
  if(length(parameters)==1) {
    f <- noquote(paste(1/diff(range(data)), "+0*exp(x)", sep=""))
  }else {
    
    f <- ToStringRe_MTE(parameters, data, 1/parExp[1]) 
  }
  
  f <- list(Function = f, Subclass = "mte",
            Domain = (fx$Domain)*sd(data)+mean(data))
  f <- new_mte(f)
  return(f)
}

#' @rdname rescaledFunctions
#' @export
ToStringRe_MTE <- function(parameters, data, num = 1)
{
  mu <- mean(data); sde=sd(data)
  str <- parameters[1]/sde
  sign <- parameters; sign[sign<0]=""; sign[sign>=0] <- "+" 
  if(mu<0) smean <- "+" else smean <- "-"
  
  if((length(parameters)-1)>0) {
    for(i in 2: length(parameters)){
      if(i<4){
        if(i%%2==1) p=-1 else p=1
        if(smean=="-") m <- -mu else m <- mu
        str <- paste(str,sign[i], (parameters[i]/sde)*exp(p*m/(num*sde)), "*exp(",ifelse(i%%2==1, "-", ""),1/(num*sde),"*x",")", sep="")
      }else{
        if(i%%2==1) p <- -(i%/%2) else p <- i%/%2
        if(smean=="-") m <- -mu else m <- mu
        str <- paste(str,sign[i], (parameters[i]/sde)*exp(p*m/(num*sde)), "*exp(",p/(num*sde),"*x)", sep="")
      }    
    }
  }
  else {
    str <- paste(str,"+0*exp(1/num*x)", sep="")
  }
  return(noquote(str))
}

#' @rdname rescaledFunctions
#' @export
meanMOP=function(fx)
{
  fx <- noquote(as.character(fx))
  t <- strsplit(fx, split="*(", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]][2]
  t1 <- strsplit(t, split=")", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]][1]
  
  mu <- strsplit(t1, split="x", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]][2]  
  mu <- as.numeric(mu)
  if(mu<0) mu <- abs(mu) else mu <- as.numeric(paste("-",mu, sep=""))
  return(mu)
}
#' @param object and object of class motbf_fit that has been learned with scaled data
#' @param POTENTIAL_TYPE the potential of the model (either "MOP" or "MTE")
#' @param data the original data set (non-scaled)
#' @rdname rescaledFunctions
#' @export
rescale_motbf_fit <- function(object, POTENTIAL_TYPE, data){
  # v_mean = sapply(data, mean)
  # v_sd = sapply(data, sd)
  v_mean = sapply(data, function(x){ifelse(is.numeric(x), mean(x, na.rm = TRUE),NA)})
  v_sd = sapply(data, function(x){ifelse(is.numeric(x), sd(x, na.rm = TRUE),NA)})
  p = length(object)
  
  # Re-scale 'functions' data.frame: variable domain and mop (or mte)
  for(j in 1:p){
    node = object[[j]]
    cpd = node$functions
    N = length(names(cpd))
    vars = names(cpd)[-N]; vars
    
    levels.node = levels(node); levels.node
    
    # re-scale domain of variables in CDP data.frame and node levels
    for(i in 1:(N-1)){
      var = vars[i]; var
      if(object[[var]]$type!="Continuous"){
        next
      }
      s = v_sd[var]; s
      m = v_mean[var]; m
      cpd[[var]] = lapply(cpd[[var]], function(x){round(x*s+m, 15)})
      levels.node[[var]] = round(levels.node[[var]]*s+m, 15)
    }
    attr(node,'levels') = levels.node
 # levels(node) = levels.node
    
    # re-scale functions of node_j
    if(object[[var]]$type=="Continuous"){
      s = v_sd[node$node]; s
      m = v_mean[node$node]; m
      
      for(i in 1:nrow(cpd)){
        attr(cpd[[N]][[i]], 'mean') = m
        attr(cpd[[N]][[i]], 'sd') = s
        f = cpd[[N]][[i]];f
        
        if(POTENTIAL_TYPE == 'MOP'){
          fx = rescaledMOP(f);fx
        }else{
          x = data[[node$node]]
          # browser()
          f$Function = gsub(paste0('*', node$node), '*x', f$Function, fixed = TRUE)
          fx = rescaledMTE(f,x)
          fx$Function = gsub('*x', paste0('*', node$node),fx$Function, fixed = TRUE)
        }
        
        cpd[[N]][[i]] = fx
      }
    }
    
    node$functions = cpd
    object[[j]] = node
  }
  
  # Re-scale levels of object
  lev = levels(object)
  vars = names(lev)
  
  for(i in 1:length(vars)){
    var = vars[i]
    if(object[[var]]$type!="Continuous"){
      next
    }
    s = v_sd[var]; s
    m = v_mean[var]; m
    lev[[var]] = round(lev[[var]]*s+m, 15)
    lev
  }
  # browser()
  attr(object, 'levels') <- lev
  # levels(object) <- lev
  return(object)
  ## End
}

#' Scale data
#' 
#' Standardize the numeric columns of a data.frame.
#' 
#' @param data data.frame containing the variables to be standardized.
#' @param v_mean vector of means of the data.frame. The value for the non-numeric columns should be NA. If not given, it is computed from data.
#' @param v_sd vector of standard deviations of the data.frame. The value for the non-numeric columns should be NA. If not given, it is computed from data.
#' @return A data.frame with scaled data. Non-numeric columns are also returned in their original position.
rescale_data= function(data, v_mean = NULL, v_sd = NULL){
  # sin recursividad
  if(is.null(v_mean)){
    v_mean = sapply(data, function(x){ifelse(is.numeric(x), mean(x, na.rm = TRUE),NA)})
  }
  if(is.null(v_sd)){
    v_sd = sapply(data, function(x){ifelse(is.numeric(x), sd(x, na.rm = TRUE),NA)})
  }
  
  n = ncol(data)
  
  res = lapply(1:n, function(i){
    x = data[[i]]
    m = v_mean[i]
    s = v_sd[i]
    if(!is.numeric(x)){
      res = x
    }else{
      res = (x-m)/s

    }
  })
  res = data.frame(res)

  colnames(res) = colnames(data)
  return(res)
}

