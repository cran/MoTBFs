#' Integrating MoTBFs
#' 
#' Compute the integral of a one-dimensional mixture of truncated basis function 
#' over a bounded or unbounded interval.
#' This function is deprecated; use \code{"integrate.motbf"} instead.
#' 
#' @param fx An object of class \code{"motbf"}.
#' @param min The lower integration limit. By default it is NULL.
#' @param max The upper integration limit. By default it is NULL.
#' @details If the limits of the interval, min and max are NULL, then the output is
#' the expression of the indefinite integral. If only 'min' contains a numeric value,
#' then the expression of the integral is evaluated at this point.
#' @return \code{integralMoTBF()} returns either the indefinite integral of the MoTBF 
#' function, which is also an object of class \code{"motbf"}, or the definite integral, 
#' wich is a \code{"numeric"} value.
#' @seealso \link{univMoTBF}, \link{integralMOP} and \link{integralMTE}
#' @export
#' 
integralMoTBF <- function(fx, min=NULL, max=NULL){
  .Deprecated("integrate.motbf")
  if(!is.motbf(fx)) stop("fx is not an 'motbf' function.")
  if(is.null(min)&&is.null(max)){
    if(is.mop(fx)) return(integralMOP(fx))
    if(is.mte(fx)) return(integralMTE(fx))
  } else if(!is.null(min)&&is.null(max)){
    if(is.mop(fx)) return(as.function(integralMOP(fx))(min))
    if(is.mte(fx)) return(as.function(integralMTE(fx))(min))
  } else{
    if(is.mop(fx)) return(as.function(integralMOP(fx))(max) - as.function(integralMOP(fx))(min))
    if(is.mte(fx)) return(as.function(integralMTE(fx))(max) - as.function(integralMTE(fx))(min))
  } 
}



# Deprecated functions

#' Integration of MOPs
#' 
#' Method to calculate the non-defined integral of an \code{"motbf"} object of \code{'mop'} subclass. 
#' This function is deprecated; use \code{"integrate.motbf"} instead.
#' 
#' @param fx An \code{"motbf"} object of subclass \code{'mop'}.
#' @return The non-defined integral of the function.
#' @seealso \link{univMoTBF} for learning and \link{integralMoTBF} 
#' for a more complete function to get defined and non-defined integrals
#' of class \code{"motbf"}.
#' @export
#' 
integralMOP <- function(fx){
  .Deprecated("integrate.motbf")
  if(!is.motbf(fx)) stop("fx is not an 'motbf' function.")
  if(is.motbf(fx)&&!is.mop(fx)) stop("fx is an 'motbf' function but not 'mop' subclass.")
  
  #options(warn=-1)
  
  suppressWarnings({
    parameters <- coeffMOP(fx)
    mu <- tryCatch(meanMOP(fx), error = function(e) NA) 
    str <- paste(parameters[1], "*x", sep="")
    if(length(parameters)==1){
      f <- noquote(str)
      f <- motbf(f)
      return(f)
    }
    for(i in 2:length(parameters)){
      if(parameters[i]>=0) sign <- "+" else sign <- ""
      if(is.na(mu)){
        str <- paste(str, sign, (parameters[i]/i), "*x^", i, sep="")
      }else{
        if(mu>=0) signmean <- "-" else signmean <- "+"
        str <- paste(str,sign,parameters[i]/i,"*(x",signmean, mu, ")^",i, sep="")
      }
    }
    f <- noquote(str)
    f <- list(Function = f, Subclass = "mop")
    f <- motbf(f)
  })
  
  return(f)
}


#' Integrating MTEs
#' 
#' Method to calculate the non-defined integral of an \code{"motbf"} object of \code{'mte'} 
#' subclass. This function is deprecated; use \code{"integrate.motbf"} instead.
#' 
#' @param fx An \code{"motbf"} object of subclass \code{'mte'}.
#' @return The non-defined integral of the function.
#' @seealso \link{univMoTBF} for learning and \link{integralMoTBF} 
#' for a more complete function to get defined and non-defined integrals
#' of class \code{"motbf"}.
#' @export
#' 
integralMTE <- function(fx){
  .Deprecated("integrate.motbf")
  if(!is.motbf(fx)) stop("fx is not an 'motbf' function.")
  if(is.motbf(fx)&&!is.mte(fx)) stop("fx is an 'motbf' function but not 'mte' subclass")
  
  parameters <- coeffMTE(fx)
  coefExponential <- coeffExp(fx)[-1]
  str <- paste(parameters[1], "*x", sep="")
  if((length(parameters)-1)>0) {
    for(i in 2:length(parameters)){
      if((parameters[i]*(1/coefExponential[i-1]))>=0) sign <- "+" else sign <- ""
      str <- paste(str,sign, parameters[i]*(1/coefExponential[i-1]), "*exp(",coefExponential[i-1],"*x)", sep="")
    }
  }else {
    ## An MTE constant
    str  <-  paste(str,"+0*exp(", 1/coefExponential[1], "*x)", sep="")
  }
  f <- noquote(str)
  f <- list(Function = f, Subclass = "mte")
  f <- motbf(f)
  return(f)
}

#' Integration with MoTBFs
#' 
#' Integrate a \code{"jointmotbf"} object over an non defined domain. It is able to
#' get the integral of a joint function over a set of variables or over all
#' the variables in the function.
#' This function is deprecated; use \code{"integrate.motbf"} instead.
#' 
#' @param P A \code{"jointmotbf"} object.
#' @param var A \code{"character"} vector containing the name of the variables that will be integrated out.
#' Instead of the names, the position of the variables can be given.
#' By default it's \code{NULL} then all the variables are integrated out.
#' @return A multiintegral of a joint function of class \code{"jointmotbf"}.
#' @export
#' 
integralJointMoTBF <- function(P, var=NULL){
  .Deprecated("integrate.motbf")
  suppressWarnings({
    Letters <- c(letters[24:26], letters[1:23])
    if(is.numeric(var)) var <- Letters[var]
    if(is.null(var)) var <- nVariables(P)
    nVar <- nVariables(P)
    for(h in 1:length(var)){
      string <- as.character(P)
      parameters <- coef(P)
      f1 <- substr(string, 1, 1)
      t <- strsplit(string, split="-", fixed = TRUE)[[1]]
      for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")##le volvemos a añadir el simbolo negativo
      if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
      if(f1 == "-"){
        t <- t[-1]
      }
      
      t2 <- c()
      for(i in 1:(length(t))){
        t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
        t2 <- c(t2,t1)
      }
      
      pos <- grep("e", t2)
      if(length(pos)!=0){
        for(i in pos){
          t2[i] <- paste(t2[i],t2[i+1], sep="")
          t2[i+1] <- NA
        }
      }
      str <- t2[!is.na(t2)]
      
      s <- strsplit(str, split=as.character(parameters), fixed=TRUE)
      for (i in 1:length(s)){
        if(is.na(s[[i]][2])) s[[i]][2]= paste("*", var[h], sep="")
        else{
          t <- strsplit(s[[i]][2], split=NULL)[[1]]
          if(!any(t%in%var[h])){
            pVar <- which(nVar%in%var[h])
            pVars <- which(nVar%in%t)
            if(any(pVar<pVars)){
              d <- which(nVar[pVars[which(pVar<pVars)[1]]]==t)
              t <- append(t, paste("*", var[h], sep=""), d-2)
              f <- c()
              for(j in 1:length(t))  f <- paste(f, t[j], sep="")
              s[[i]][2] <- f
            }else{
              s[[i]][2] <- paste(s[[i]][2], "*", var[h], sep="")
            }
          }else{
            p <- which(t%in%var[h])
            if(is.na(as.numeric(t[p+2]))){
              t <- append(t, paste("^" , 2, sep=""),after=p)
              parameters[i] <- parameters[i]/2
            }else{
              t[p+2]=as.numeric(t[p+2])+1
              parameters[i] <- parameters[i]/as.numeric(t[p+2])
            } 
            f <- c()
            for(j in 1:length(t))  f <- paste(f, t[j], sep="")
            s[[i]][2] <- f
          }
          
        }
      }
      if(any(parameters[2:length(parameters)]>=0)) parameters[parameters[2:length(parameters)]>=0]
      parameters[2:length(parameters)][parameters[2:length(parameters)]>=0] = 
        paste("+",parameters[2:length(parameters)][parameters[2:length(parameters)]>=0], sep="")
      
      s <- unlist(s)
      s[which(s%in%"")]=parameters
      str <- c()
      for(i in 1:length(s)) str <- paste(str, s[i], sep="")
      
      str <- noquote(str)
      if(length(nVariables(str))>1) P <- jointmotbf(str)
      else P <- motbf(str)
    }
    
    return(P)
  })
}






#' Evaluation of joint MoTBFs
#' 
#' Evaluates a \code{"jointmotbf"} object at a specific point. 
#' This function is deprecated; use \code{"eval.mop"} instead.
#' 
#' @param P A \code{"jointmotbf"} object.
#' @param values A list with the name of the variables equal to the values to be evaluated.
#' @return If all the variables in the equation are evaluated then a \code{"numeric"} value
#' is returned. Otherwise, an \code{"motbf"} object or a \code{"jointmotbf"} object is returned.
#' @export
#' 
evalJointFunction <- function(P, values){
  .Deprecated("eval.mop")
  
  nVar <- nVariables(P)
  if(is.motbf(P)) return(as.numeric(as.function(P)(values[[1]])))
  string <- as.character(P)
  parameters <- coef(P)
  f1 <- substr(string, 1, 1)
  t <- strsplit(string, split="-", fixed = TRUE)[[1]]
  for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")##le volvemos a añadir el simbolo negativo
  if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
  if(f1 == "-"){
    t <- t[-1]
  }
  
  t2 <- c()
  for(i in 1:(length(t))){
    t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  
  pos <- grep("e", t2)
  if(length(pos)!=0){
    for(i in pos){
      t2[i] <- paste(t2[i],t2[i+1], sep="")
      t2[i+1] <- NA
    }
  }
  str <- t2[!is.na(t2)]
  
  s <- strsplit(str, split=as.character(parameters), fixed=TRUE)
  var <- names(values)
  for(i in 1:length(s)){
    for(j in 1:length(var)){
      if(length(grep(var[j],s[[i]]))==0){
        if(j==1){s[[i]][1]=parameters[i]; if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")}
      }else{
        if(j==1){s[[i]][1]=parameters[i]; if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")}
        if(!is.na(s[[i]][2])){
          s1=strsplit(s[[i]][2], split=NULL)[[1]]
          if(any(s1%in%var[j])){
            if(!is.na(s1[which(s1==var[j])+1])){
              if(s1[which(s1==var[j])+1]=="^"){
                exp=as.numeric(s1[which(s1==var[j])+2])
                s[[i]][1]=as.numeric(s[[i]][1])*(values[[j]]^exp)
                if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
                s1[which(s1==var[j])+1]=""
                s1[which(s1==var[j])-1]=""
                s1[which(s1==var[j])+2]=""
                s1[s1%in%var[j]]=""
                s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
                s[[i]][2]=s2
              }else{
                s[[i]][1]=as.numeric(s[[i]][1])*values[[j]]
                if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
                s1[which(s1%in%var[j])-1]=""
                s1[which(s1%in%var[j])]=""
                s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
                s[[i]][2]=s2
              }
            }else{
              s[[i]][1]=as.numeric(s[[i]][1])*values[[j]]
              if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
              s1[which(s1%in%var[j])-1]=""
              s1[which(s1%in%var[j])]=""
              s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
              s[[i]][2]=s2
            }
          } 
        }
      }
    }
  }
  
  if(is.na(s[[1]][2])) s[[1]][2]=""
  exp <- sapply(1:length(s), function(i) s[[i]][2])
  exp <- unique(exp)  
  s <- unlist(s)
  param <- c()
  for(i in 1:length(exp)) param <- c(param,sum(as.numeric(s[which(s==exp[i])-1])))
  sign <- ""
  if(length(param)>1) for(i in 2:length(param)) if(param[i]>0) sign <- c(sign,"+") else sign <- c(sign,"")
  str <- c()
  for(i in 1:length(exp)) str <- paste(str, sign[i], param[i], exp[i], sep="")
  
  if(all(nVar%in%names(values))){
    message("The method 'as.function()' can be used. \n")
    return(eval(parse(text=str)))
  } else{
    f <- noquote(str)
    l <- length(nVariables(f))
    if(l==1) f <- motbf(f)
    if(l>1) f <- jointmotbf(f)
    return(f)
  }  
}





#' Joint MoTBF density learning
#' 
#' Deprecated functions; use \code{"jointmotbf.fit"} instead.
#' 
#' Two functions for learning joint MoTBFs. 
#' The first one, \code{parametersJointMoTBF()},
#' gets the parameters by solving a quadratic optimization problem, minimizing
#' the mean squared error between the empirical joint CDF and the estimated CDF.
#' The density is obtained as the derivative od the estimated CDF.
#' The second one, \code{jointMoTBF()}, fixes the equation of the joint function using 
#' the previously learned parameters and converting this \code{"character"} string into an 
#' object of class \code{"jointmotbf"}.
#' 
#' @name jointmotbf.learning
#' @rdname jointmotbf.learning
#' @param X A dataset of class \code{"data.frame"}.
#' @param ranges A \code{"numeric"} matrix containing the range of the varibles used to fit the function, 
#' where each column corresponds to a variable. If not specified, the range of each variable is computed from the data.
#' @param dimensions A \code{"numeric"} vector containing the number of parameters of each varible.
#' @param object A list with the output of the function \code{parametersJointMoTBF()}.
#' @return
#' \code{parametersJointMoTBF()} returns a list with the following elements: 
#' \bold{Parameters}, which contains the computed coefficients of the resulting function; 
#' \bold{Dimension}, which is a \code{"numeric"} vector containing the number 
#' of coefficients used for each variable; 
#' \bold{Range} contains a \code{"numeric"} matrix with the domain of each variable, by columns;
#' \bold{Iterations} contains the number of iterations needed to solve the problem;
#' \bold{Time} contains the execution time.
#' 
#' \code{jointMoTBF()} returns an object of class \code{"jointmotbf"}, which is a list whose only visible element
#' is the analytical expression of the learned density. It also contains the other aforementioned elements, 
#' which can be retrieved using \code{attributes()}
#' @importFrom Matrix nearPD
#' @export
parametersJointMoTBF=function(X, ranges=NULL, dimensions=NULL){
  .Deprecated("jointmotbf.fit")
  
  tm <- Sys.time()
  
  if(is.null(ranges)) ranges <- sapply(X, range)
  if(is.null(dimensions)){
    Fx <- lapply(X, univMoTBF, "MOP")
    ll <- list(); for(i in 1:length(Fx)) ll[[length(ll)+1]] <- length(coef(Fx[[i]]))
  } else {
    ll <- list(); for(i in 1:length(dimensions)) ll[[length(ll)+1]] <- dimensions[i]
  }
  
  dim <- c(); pos <- which(c(length(ll[[1]]), length(ll[[2]]))==min(length(ll[[1]]), length(ll[[2]])))
  if(pos[1]==1) {pos1 <- 1; pos2 <- 2}
  if(pos[1]!=1) {pos1 <- 2; pos2 <- 1}
  for(i in 1:length(ll[[pos1]])) for(j in 1:length(ll[[pos2]])) dim <- c(dim,list(c(ll[[pos1]][i],ll[[pos2]][j])))
  if(pos1==2) for(i in 1:length(dim)) dim[[i]] <- dim[[i]][length(dim[[i]]):1]
  
  size <- prod(sapply(1:length(ll), function(i) length(ll[[i]])))
  si <- lapply(1:ncol(X),function(i) rep(ll[[i]][1], size))
  dim <- c(); for(i in 1:length(si)) dim <- cbind(dim,unlist(si[[i]]))
  dim <- lapply(1:nrow(dim), function(i) dim[i,])
  
  ## Create the grid: depending on the n.record and the n.variables
  #npointsgrid <- 100
  # npointsgrid <- max(ceiling((2/(ncol(X)))^2*100), 10)
  # npointsgrid <- 10
  fitPoints <- max(ceiling((2/(ncol(X)))^2*100), 10)
  
  eg <- list()
  for(i in 1: ncol(X)){
    eg[[i]] <- seq(ranges[1,i], ranges[2,i], length.out = fitPoints)
    
  }
  x <- expand.grid(eg)
  
  ## Cumulative densities
  y <- jointCDF(X,x)
  
  P <- c()
  for(s in 1:length(dim)){
    Xt <- ranges; pos <- coefExpJointCDF(dim[[s]])
    n <- nrow(x); nterms <- length(pos)
    
    # x1 <- c()
    # for(j in 1:ncol(X)){
    #   xx=rep(1, n)
    #   for(i in 1:nterms){
    #     xi <- cbind((x[,j]^pos[[i]][j]))
    #     xx <- cbind(xx,xi)
    #   }
    #   x1[[length(x1)+1]] <- xx
    # }
    
    x1 <- list()
    for(j in 1:ncol(X)){
      pos2 <- c(0,unlist(lapply(pos, function(x){x[j]})))
      x1[[j]] <- outer(x[,j], pos2, FUN = "^")
    }
    xx <-Reduce("*", x1)
    
    ## Dmat
    XX <- t(xx)%*%xx 
    
    ## dvec
    Xy <- t(xx)%*%y 
    
    ## Constrains
    xNew <- list()
    constraints = round(60000^(1/ncol(Xt)))
    
    for( j in 1:ncol(Xt)){
      # xnew <- seq(min(Xt[,j]), max(Xt[,j]), length=round(60000^(1/ncol(Xt))))
      xnew <- seq(min(Xt[,j]), max(Xt[,j]), length=constraints)
      if(xnew[length(xnew)]!=max(Xt[,j])) xnew <- c(xnew,max(Xt[,j]))
      xNew[[length(xNew)+1]] <- xnew
    }
    xN <- expand.grid(xNew)
    
    ma22 <- c()
    for(i in 1:ncol(X)){
      d <- sapply(1:nterms, function(j) pos[[j]][i]*(xN[,i]^(pos[[j]][i]-1)))
      tf <- sapply(1:length(pos), function(j) pos[[j]][i]==0)
      if(any(tf)) d[,which(tf)] <- 0
      d <- cbind(rep(0, nrow(xN)),d)
      ma22[[length(ma22)+1]] <- t(d)
    }
    ma2 <- Reduce("*", ma22)
    
    ma3 <- c()
    for(j in 1:ncol(Xt)){
      ma <- 0
      for(i in 1:nterms) ma <- rbind(ma,(max(Xt[,j])^pos[[i]][j]) - (min(Xt[,j])^pos[[i]][j]))
      ma3 <- cbind(ma3,ma)
    }
    for(i in 1:(ncol(X)-1)) ma3[,i+1] <- ma3[,i] * ma3[,i+1]
    ma33 <- cbind(ma3[,i+1])
    
    ## Amat
    AA <- cbind(ma33, ma2)
    
    ## bvec
    B <- c(1,rep(1.0E-5,ncol(ma2)))
    
    tr <- tryCatch(solve.QP(XX, Xy, AA, B, meq=1), error = function(e) NULL)
    if(is.null(tr)==TRUE){
      message("matrix D in quadratic function has been approximated to the nearest positive definite!")
      tr <- solve.QP(nearPD(XX)$mat, Xy, AA, B, meq=1)
      # return(NULL)
    }   
    
    finaltm <- Sys.time() - tm
    soluc <- tr 
    parameters <- soluc$solution
    
    # if(is.null(tr)==TRUE){
    #   return(NULL)
    # }else{    
    #   soluc <- tr 
    #   parameters <- soluc$solution
    # }
    
    ## Parameters PDF
    a <- sapply(pos, prod)
    param <- parameters[-1]*a
    parameters <- param[which(param!=0)]
    
    P <- list(Parameters = parameters, 
              Dimension = dim[[s]], Range = Xt, 
              Iterations = tr$iterations[1], 
              Time=finaltm)
  }
  
  return(P)
}

#' @rdname jointmotbf.learning
#' @export
jointMoTBF <- function(object){
  .Deprecated("jointmotbf.fit")
  
  parameters <- object$Parameters
  dimensions <- object$Dimension
  v <- letters[24:(24+length(dimensions)-1)]
  v[is.na(v)] <- letters[1:length(v[is.na(v)])]
  if(length(parameters)==1) {
    s <- paste(parameters[1], "+0", sep="")
    for(i in 1:length(v)) s <- paste(s, "*", v[i], sep="")
    P <- list(Function = s, Domain = object$Range)
    P <- jointmotbf(P)
    return(P)
  }
  pos <- coefExpJointCDF(dimensions)
  t <- list()
  for(i in 1:length(pos)) if(prod(pos[[i]])!=0) 
    t[[length(t)+1]] <- pos[[i]] else next
  pos <- t
  
  if(any(dimensions==1)) posvardim1 <- which(dimensions==1)
  
  if(parameters[1]==0) s <- c() else s <- parameters[1]
  for(i in 2:length(parameters)){
    if(is.null(s)){
      s <- paste(s, parameters[i],sep="")
    } else{
      s <- paste(s, ifelse(parameters[i]>=0, "+", ""), parameters[i],sep="")
      if(i != length(parameters)) for(p in 1:length(v)) s <- paste(s, ifelse((pos[[i]][p]-1)==0, "", paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "", paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
      if(i == length(parameters)){
        for(p in 1:length(v)) s <- paste(s, ifelse((pos[[i]][p]-1)==0, ifelse(any(p==posvardim1),
                                                                              paste("*",v[p],"^",pos[[i]][p]-1, sep=""),""), 
                                                   paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "",
                                                                          paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
      }
    }
  }
  P <- noquote(s)
  P <- list(Function = P, Domain = object$Range,
            Iterations = object$Iterations, Time = object$Time)
  P <- jointmotbf(P)
  return(P)
}




#' Dimension of MoTBFs
#' 
#' Get the dimension of \code{"motbf"} and \code{"jointmotbf"} densities.
#' This function is deprecated; use \code{"getMotbfDim"} instead.
#' 
#' @param P An object of class \code{"motbf"} and subclass 'mop' or \code{"jointmotbf"}.
#' @return Dimension of the function.
#' @seealso \link{univMoTBF} and \link{jointMoTBF}
#' @export
#'
dimensionFunction <- function(P){
  .Deprecated("getMotbfDim")
  suppressWarnings({
    string <- as.character(P)
    nVar <- nVariables(P)
    
    if(length(nVar)==1){
      str <- strsplit(string, split="^", fixed=TRUE)[[1]]
      return(as.numeric(str[length(str)])+1)
    }
    dimensions <- c()
    for(i in 1:length(nVar)){
      elements <- strsplit(string, split=NULL)[[1]]
      pos <- grep(nVar[i], elements)
      nele <- elements[(pos[length(pos)]+1):(pos[length(pos)]+3)]
      if(is.na(as.numeric(nele[3]))) nele <- nele[-3]
      if(length(which(nele%in%"^"))!=0){
        #nele[1]=="^")
        pos <- as.numeric(nele[c(2,3)])
        if(all(!is.na(pos))){
          num <- c()
          for(h in 2:length(nele)) num <- paste(num,nele[h], sep="")
          dimensions=c(dimensions, as.numeric(num))
        }else dimensions=c(dimensions, as.numeric(nele[2]))
      } else dimensions=c(dimensions, 1)
    }
    return(dimensions+1)
  })
}




#' Number of Variables in a Joint Function
#' 
#' Compute the number of variables which are in a \code{jointmotbf} object.
#' This function is deprecated; use \code{"getMotbfVar"} instead.
#' 
#' @param P An \code{"motbf"} object or a \code{"jointmotbf"} object.
#' @return A \code{"character"} vector with the names of the variables in the function.
#' @export
#' 
nVariables=function(P){
  .Deprecated("getMotbfVar")
  suppressWarnings({
    Letters <- letters[c(24:26,1:23)]
    string <- as.character(P)
    t <- strsplit(string, split="-", fixed = TRUE)[[1]]
    t2 <- c()
    for(i in 1:length(t)){
      t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
      t2 <- c(t2,t1)
    }
    t <- unlist(sapply(1:length(t2), function(i) strsplit(t2[i], split="*", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]))
    
    variables <- t[is.na(as.numeric(t))]
    t <- unlist(strsplit(variables, split="^", fixed = TRUE, perl = FALSE, useBytes = FALSE))
    t <- unlist(strsplit(t, split="(", fixed = TRUE, perl = FALSE, useBytes = FALSE))
    variables <- t[is.na(as.numeric(t))]
    variables <- unique(variables)
    variables <- Letters[Letters%in%variables]
  })
  return(variables)
}



#' Marginalization of MoTBFs
#'
#' Computes the marginal densities from a \code{"jointmotbf"}
#' object. 
#' This function is deprecated; use \code{"marginal.jointmotbf"} instead.
#'
#'@param P An object of class \code{"jointmotbf"}, i.e., the joint density function.
#'@param var The \code{"numeric"} position or the \code{"character"} name of the marginal variable.
#'@return The marginal of a \code{"jointmotbf"} function. The result is an object of class \code{"motbf"}.
#'@seealso \link{jointMoTBF} and \link{evalJointFunction}
#'@export
#' 
marginalJointMoTBF <- function(P, var){
  
  .Deprecated("marginal.jointmotbf")
  
  suppressWarnings({
    Letters <- c(letters[24:26], letters[1:23])
    if(is.numeric(var)) var <- Letters[var] else var <- var
    varN <- nVariables(P); noVar=varN[!varN%in%var]
    l <- length(noVar)
    for(j in 1:l){
      varN <- nVariables(P)
      posnoVar <- which(varN%in%noVar[j])
      posvar <- which(varN%in%var)
      
      min <- lapply(posnoVar, function(i) P$Domain[1,i])
      names(min) <- noVar[j]
      max <- lapply(posnoVar, function(i) P$Domain[2,i])
      names(max) <- noVar[j]
      
      iP <- integralJointMoTBF(P, noVar[j])
      minF <- evalJointFunction(iP, min)
      maxF <- evalJointFunction(iP, max)
      if(is.numeric(minF)&&is.numeric(maxF)){
        P <- list(Function=(maxF-minF), Subclass="mop", Domain= P$Domain[,which(!varN%in%noVar[j])])
        P <- motbf(P)
        return(P)
      }
      
      ## Min part
      nVar <- nVariables(minF)
      string <- as.character(minF)
      f1 <- substr(string, 1, 1)
      t <- strsplit(string, split="-", fixed = TRUE)[[1]]
      for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")
      if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
      if(t[1]=="-") t=t[-1]
      t2 <- c()
      for(i in 1:(length(t))){
        t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
        t2 <- c(t2,t1)
      }
      pos <- grep("e", t2)
      if(length(pos)!=0){
        for(i in pos){
          t2[i] <- paste(t2[i],t2[i+1], sep="")
          t2[i+1] <- NA
        }
      }
      str <- t2[!is.na(t2)]
      coefMIN <- coef(minF)
      
      sMIN <- unlist(strsplit(str, split=coefMIN, fixed=TRUE))
      varseq=(unique(sMIN))[-1]
      sMIN[which(sMIN=="")]=coefMIN
      paramMIN <-c()
      for(i in 1:length(varseq)){
        pos <- which(sMIN%in%varseq[i])
        paramMIN <- c(paramMIN, sum(as.numeric(sMIN[pos-1])))
        sMIN <- sMIN[-c(pos, pos-1)]
      }
      paramMIN <- c(sum(as.numeric(sMIN)), paramMIN)
      
      ## Max part
      nVar <- nVariables(maxF)
      string <- as.character(maxF)
      f1 <- substr(string, 1, 1)
      t <- strsplit(string, split="-", fixed = TRUE)[[1]]
      for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")
      if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
      if(t[1]=="-") t=t[-1]
      t2 <- c()
      for(i in 1:(length(t))){
        t1 <- strsplit(t[i], split="+", fixed = TRUE, perl = FALSE, useBytes = FALSE)[[1]]
        t2 <- c(t2,t1)
      }
      pos <- grep("e", t2)
      if(length(pos)!=0){
        for(i in pos){
          t2[i] <- paste(t2[i],t2[i+1], sep="")
          t2[i+1] <- NA
        }
      }
      str <- t2[!is.na(t2)]
      coefMAX <- coef(maxF)
      sMAX <- unlist(strsplit(str, split=coefMAX, fixed=TRUE))
      varseq=(unique(sMAX))[-1]
      sMAX[which(sMAX=="")]=coefMAX
      paramMAX <-c()
      for(i in 1:length(varseq)){
        pos <- which(sMAX%in%varseq[i])
        paramMAX <- c(paramMAX, sum(as.numeric(sMAX[pos-1])))
        sMAX <- sMAX[-c(pos, pos-1)]
      }
      paramMAX <- c(sum(as.numeric(sMAX)), paramMAX)
      
      ## New parameters
      newparam <- paramMAX - paramMIN
      
      sign=""
      for(v in 2:length(newparam)) if(newparam[v]>0) sign=c(sign,"+") else sign=c(sign,"")
      st <- newparam[1]
      for(v in 2:length(newparam)) st <- paste(st, sign[v],newparam[v], varseq[v-1], sep="")
      if(length(nVariables(st))>1){
        P <- list(Function=noquote(st), Domain= P$Domain[,which(!varN%in%noVar[j])])
        P <- jointmotbf(P)
      }else{
        P <- list(Function=noquote(st), Subclass="mop", Domain= P$Domain[,which(!varN%in%noVar[j])])
        P <- motbf(P)
        if(length(unique(coeffPol(P)))==1){
          st <- paste(sum(coef(P)), ifelse(unique(coeffPol(P))==0, "",paste("*", var, ifelse(unique(coeffPol(P))==1, "", paste("^", unique(coeffPol(P)), sep="")), sep="")), sep="")
          P <- list(Function=noquote(st), Subclass="mop", Domain= P$Domain)
          P <- motbf(P)
        }
      }
    }
    return(P)
  })
}