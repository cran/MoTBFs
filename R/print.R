#' Print object of class motbf
#' \code{print} method for class \code{"motbf"}.
#' @param x An object of class \code{"motbf"}, \code{"motbf_fit"}, \code{"univmotbf"}, \code{"piecewisemop"}, \code{"jointmotbf"}, or \code{"motbf.fit.cv"}.
#' @param ... optional arguments passed to print for other classes created in the MoTBFs package. 
#' Currently, no optional arguments are supported.
#' @details The following classes are created in the MoTBFs package:
#' \describe{
#' \item{motbf}{the generic class common to all objects created in the package}
#' \item{univmotbf}{the class corresponding to MTE or MOP univariate distributions. An object of class \code{univmotbf} is the output of function \code{univMoTBF()}.}
#' \item{piecewisemop}{the class corresponding to MOP univariate distributions defined by multiple sub-functions with different domain. When calling \code{variableElimination()}, the output might be of this class.}
#' \item{jointmotbf}{the class corresponding to joint distributions. An object of class \code{jointmotbf} is the output of function \code{jointMOP()}}
#' \item{motbf_fit}{the class corresponding to fully fitted Bayesian network models (either discrete, continuous or hybrid). An object of class \code{motbf_fit} is the output of function \code{motbf.fit()}.}
#' \item{motbf.fit.cv}{the class corresponding to k-fold cross validation results. An object of class \code{motbf.fit.cv} is the output of function \code{motbf.cv()}.}
#' \item{motbf.fit.node}{the class corresponding to a single node of a Bayesian Network.}
#' }
#' @rdname print.motbf
#' @exportS3Method base::print motbf
#' @examples
#' \donttest{
#' ## Dataset Ecoli
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' P <- motbf.fit(graph = dag, data = data, numIntervals = 3, POTENTIAL_TYPE = "MOP",
#' maxParam = 15)
#' P
#' }
#' 
print.motbf <- function(x, ...){ 
  NextMethod("print")
}


#' @exportS3Method base::print motbf_fit
#' @rdname print.motbf
print.motbf_fit <- function(x, ...){
  if(!("motbf_fit"%in%class(x))){
    stop("The object is not of class motbf_fit")
  }
  
  for(i in 1:length(x)){
    node = x[[i]]
    fx = node$functions
    n = ncol(fx)
    
    if(is.null(node$parents)){
      cat("Potential(", node$node, ")\n", sep="")
      if(node$type == 'Discrete'){
        fx = collapseDiscreteCPD(node)
        dimnames = list(c(),unlist(fx[[1]]))
        print(as.data.frame(matrix(fx[[1,n]], dimnames = dimnames, nrow = 1)), row.names = FALSE)
        
      }else{
        print(fx[[1,ncol(fx)]])
      }
      cat("\n")
    }else{
      cat("Potential(", node$node," | ",paste0(node$parents, collapse = ", "),")\n", sep="")
      # printConditionalBN(node)
      print.motbf.fit.node(node)
      cat("\n")
    } 
    
  }
}

#' Print a single node of a BN. 
#' This function is called by print.motbf_fit, but not exported.
#' @rdname print.motbf
#' @exportS3Method base::print motbf.fit.node
print.motbf.fit.node <- function(x,...){
  node = x
  fx = node$functions
  n = ncol(fx)
  pa = node$parents
  X = node$node
  
  if(is.null(pa)){ # node with no parents
    cat("Potential(", node$node, ")\n", sep="")
    if(node$type == 'Discrete'){
      fx = collapseDiscreteCPD(node)
      dimnames = list(c(),unlist(fx[[1]]))
      print(as.data.frame(matrix(fx[[1,n]], dimnames = dimnames, nrow = 1)), row.names = FALSE)
      
    }else{# continuous node with no parents
      print(fx[[1,ncol(fx)]])
    }
    # End of function for this case (node with no parents)
    
  }else{# node with parents
    # states = lapply(fx[,c(pa,X)], unique)
    states = attr(node, 'levels')
    p_dist = fx[,-c(n,(n-1)), drop = FALSE]
    
    if(all(sapply(states,is.character))){# both node and parents are discrete
      pa = colnames(fx)[(n-2):1, drop = FALSE]
      
      param = fx[,n]
      
      # cptX_dimnames = lapply(fx[,c(X,pa)], unique)
      cptX_dimnames = states[c(X,pa)]
      
      cptX_dim = sapply(cptX_dimnames, length)
      
      if(nrow(p_dist)!=prod(cptX_dim)){ # This piece of code is no longer necessary. CPTs for discrete nodes are fixed
        combinations = expand.grid(states[names(cptX_dimnames)])[,-1, drop = FALSE]
        if(ncol(combinations)>1){
          combinations = combinations[,ncol(combinations):1]
        }
        
        A = apply(p_dist, 1, paste, collapse = ',')
        B = apply(combinations, 1, paste, collapse = ',')
        
        combNotObs = unique(B[which(!(B%in%A))])
        
        if(length(combNotObs)>0){
          # get uniform distribution for combinations not observed in the data
          unifDistNotObs = 1/length(which(B==combNotObs[1])) 
          # logic vector: T if combinations are observed; F if not
          cptLogic = B%in%A
          # keep CPT where combinations are observed
          cptLogic[which(cptLogic == TRUE)] = param
          # add uniform distribution where combinations are not observed
          cptLogic[!(B%in%A)] = unifDistNotObs
          
          param = cptLogic
        }
      }
      
      
      print(array(param, dim = cptX_dim, dimnames = cptX_dimnames))
      # End of function for this case (node and parents are discrete)
      
    }else{
      if(node$type == 'Discrete'){# node is discrete and at least one parent is continuous
        fx = collapseDiscreteCPD(node)
        dimnames = list(c(),unique(unlist(fx[[node$node]])))
        cpd = as.data.frame(matrix(unlist(fx[[n]]), byrow = TRUE, nrow = nrow(fx), dimnames = dimnames))
        cpd = unlist(apply(cpd, 1, list), recursive = FALSE)
        
      }else{# node is continuous
        cpd = fx[[n]]
        # cpd = lapply(cpd, function(x){
        #   unname(data.frame(x$Function))
        # })
      }
      
      # loop over cases (combinations of parents)
      for(i in 1:nrow(fx)){
        
        # loop over parents
        for(k in 1:length(pa)){
          if(is.numeric(fx[[i, pa[k]]])){
            cat("Parent:", pa[k], "  \t Range:", fx[[i, pa[k]]][1], "<", pa[k], "<", fx[[i, pa[k]]][2],'\n')
          }else{
            cat("Parent:", pa[k], "  \t Value = ", paste("\"",fx[[i, pa[k]]],"\"", sep=""),'\n')
          }
          
        }
        
        
        print(cpd[[i]], row.names = FALSE, quote = FALSE)
        # cat('\n')
      }
    }
    
  }
  
}



#' Print the results of a k-fold cross validation
#' @rdname print.motbf
#' @exportS3Method base::print motbf_fit_cv
print.motbf_fit_cv <- function(x,...){
  
  k = attr(x, 'k')
  fold.loss = sapply(x, '[[', 'loss')
  loss = attr(x, 'loss.info')
  loss.matrix = attr(x, 'loss.matrix')
  if(loss == 'logl'){
    loss = 'Log-likelihood'
  }else if(loss == 'pred'){
    if(attr(x, 'target.type') == 'Discrete'){
      if(!is.null(loss.matrix)){
        loss = 'Weighted classification accuracy'
      }else{
        loss = 'Classification accuracy'
      }
    }else{
      loss = 'Root mean squared error'
    }
  }
  if(k == 0){
    cat('Validation using training set\n', sep = '')
  }else if (k == 1){
    cat('Hold-out validation \n', sep = '')
  }else{
    cat(k, '-fold cross validation \n', sep = '')
  }
  
  cat(' Loss:', loss, '\n')
  cat(' Loss value by fold:', fold.loss, '\n')
  cat(' Average loss value:', mean(fold.loss), '\n')
}

#' @rdname print.motbf
#' @exportS3Method base::print univmotbf
print.univmotbf <- function(x,...){
  print(x[[1]])
}

#' @rdname print.motbf
#' @exportS3Method base::print piecewisemop
print.piecewisemop = function(x, ...){
  sapply(1:length(x), function(i){
    cat(noquote(paste0(' Domain: (', paste0(as.vector(x[[i]]$Domain), collapse = ', '), ')\n')))
    cat('',noquote(x[[i]]$Function))
    cat('\n')
  })
}



#' @rdname print.motbf
#' @exportS3Method base::print jointmotbf
print.jointmotbf <- function(x, ...) print(x[[1]])




## SUMMARY AND PRINT SUMMARY ####



#' Summarize an \code{"motbf"} object by describing its main features.
#' 
#' @name summary.motbf
#' @param object An object of class \code{"motbf"}.
#' @param x An object of class \code{"summary.motbf"}.
#' @param ... further arguments passed to or from other methods.
#' @return The summary of an \code{"motbf"} object. It contains a list of
#' elements with the most important information of the object.
#' @seealso \link{univMoTBF}
#' @examples
#' ## Subclass 'MOP'
#' X <- rnorm(1000)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MOP") ## or POTENTIAL_TYPE="MTE"
#' summary(P)
#' attributes(sP <- summary(P))
#' attributes(sP)
#' sP$Function
#' sP$Subclass
#' sP$Iterations
#' 
#' ## Subclass 'MTE'
#' X <- rnorm(1000)
#' P <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' summary(P)
#' attributes(sP <- summary(P))
#' attributes(sP)
#' sP$Function
#' sP$Subclass
#' sP$Iterations
#' @rdname summary.motbf
#' @exportS3Method base::summary motbf
#' 
summary.motbf <- function(object, ...){
  NextMethod("summary")
}



#' @rdname summary.motbf
#' @exportS3Method base::summary motbf_fit_cv
summary.motbf_fit_cv <- function(object, ...){
  x = object
  k = length(x)
  target.type = attr(x, 'target.type')
  loss = attr(x, 'loss.info')
  if(loss == 'pred'){
    if(target.type == 'Discrete'){
      res = confusionMatrix(x)
    }else{
      res = gof(x, method = 'pred',...)
    }
  }else{
    res = gof(x, method = 'logl')
  }
  
  return(res)
}

#' @rdname summary.motbf
#' @exportS3Method base::summary univmotbf
summary.univmotbf <- function(object, ...){
  l <- lapply(1:length(object), function(i) object[[i]])
  names(l) <- attributes(object)$names
  l[[length(l)+1]] <- class(object)
  names(l)[length(l)] <- "Class"
  l[[length(l)+1]] <- coef(object)
  names(l)[length(l)] <- "Coeff"
  class(l) <- c("summary.univmotbf", "motbf")
  l
}

#' @rdname summary.motbf
#' @exportS3Method base::print summary.univmotbf
print.summary.univmotbf <- function(x, ...)
{
  cat("\n MoTBFs FOR UNIVARIATE DISTRIBUTIONS \n")
  cat("\n Model:"); cat("\n", x$Function, "\n")
  cat("\n Class:", x$Class)
  cat("\n Subclass:", x$Subclass, "\n")
  cat("\n Coefficients:"); cat("\n", x$Coeff, "\n")
  if(!is.null(x$Domain)){ 
    cat("\n Domain:")
    cat("\n (",x$Domain[[1]], ", ",x$Domain[[2]],")", sep="","\n")
  }
  if(!is.null(x$Iterations)&&!is.null(x$Time)){
    cat("\n Number of Iterations:", x$Iterations, "\n")
    cat("\n Processing Time:", x$Time, attributes(x$Time)$units, "\n")
  }
  invisible(x)
}

#' @rdname summary.motbf
#' @exportS3Method base::summary piecewisemop
summary.piecewisemop <- function(object, ...)
{
  
  l <- lapply(1:length(object), function(i) object[[i]])
  l <- list(object)
  names(l) <- 'Functions'
  l[[length(l)+1]] <- class(object)
  names(l)[length(l)] <- "Class"
  l[[length(l)+1]] <- sapply(object, coef)
  names(l)[length(l)] <- "Coeff"
  class(l) <- c("summary.piecewisemop", "motbf")
  l
}

#' @rdname summary.motbf
#' @exportS3Method base::print summary.piecewisemop
print.summary.piecewisemop <- function(x, ...)
{
  cat("\n MoTBFs FOR UNIVARIATE DISTRIBUTIONS \n")
  cat("\n Model: \n")
  
  print(x$Functions)
  
  cat("\n Class:", x$Class)
  cat("\n Coefficients:"); cat("\n", x$Coeff, "\n")
  if(!is.null(x$Domain)){ 
    cat("\n Domain:")
    cat("\n (",x$Domain[[1]], ", ",x$Domain[[2]],")", sep="","\n")
  }
  if(!is.null(x$Iterations)&&!is.null(x$Time)){
    cat("\n Number of Iterations:", x$Iterations, "\n")
    cat("\n Processing Time:", x$Time, attributes(x$Time)$units, "\n")
  }
  invisible(x)
}


#' @rdname summary.motbf
#' @exportS3Method base::summary jointmotbf
summary.jointmotbf <- function(object, ...)
{
  l <- lapply(1:length(object), function(i) object[[i]])
  names(l) <- attributes(object)$names
  l[[length(l)+1]] <- class(object)
  names(l)[length(l)] <- "Class"
  l[[length(l)+1]] <- coef(object)
  names(l)[length(l)] <- "Coeff"
  class(l) <- c("summary.jointmotbf", "motbf")
  l
}

#' @rdname summary.motbf
#' @exportS3Method base::print summary.jointmotbf
print.summary.jointmotbf <- function(x, ...)
{
  cat("\n MoTBFs FOR MULTIVARIATE DISTRIBUTIONS \n")
  cat("\n Model:"); cat("\n", x$Function, "\n")
  cat("\n Class:", x$Class, "\n")
  cat("\n Coefficients:"); cat("\n", x$Coeff, "\n")
  if(!is.null(x$Domain)){
    varname = colnames(x$Domain)
    for(i in 1:ncol(x$Domain)){
      cat("\n Domain ",varname[i], ":", sep="")
      cat("\n (",x$Domain[1,i], ", ",x$Domain[2,i],")", sep="")
    }
    cat("\n")
  }
  if(!is.null(x$Iterations)&&!is.null(x$Time)){
    cat("\n Number of Iterations:", x$Iterations, "\n")
    cat("\n Processing Time:", x$Time, attributes(x$Time)$units, "\n")
  }
  invisible(x)
}


