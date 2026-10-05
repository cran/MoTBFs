#' Learning conditional MoTBF densities
#' 
#' Collection of functions used for learning conditional MoTBFs,
#' computing the internal BIC, selecting the parents that get 
#' the best BIC value, and other internal functions required to learn the
#' conditional densities.
#' 
#' @name conditionalmotbf.learning
#' @rdname conditionalmotbf.learning
#' 
#' @param data An object of class \code{"data.frame"}.
#' @param nameParents A \code{"character"} vector containing the names of the parent variables.
#' @param nameChild A \code{"character"} string containing the name of the child variable.
#' @param domainChild A \code{"numeric"} vector with the range of the child variable.
#' @param domainParents An object of class \code{"matrix"} with the range of the parent 
#' variables, or a \code{"numeric"} vector if there is only one parent.
#' @param numIntervals A positive integer indicating the maximum number of intervals 
#' for splitting the domain of the parent variables.
#' @param POTENTIAL_TYPE A \code{"character"} string, either \emph{MOP} or \emph{MTE}, corresponding 
#' to the type of basis function.
#' @param maxParam A positive integer which indicates the maximum number of coefficients in the 
#' function. If specified, the output is the function which gets the best BIC with, at most, 
#' this number of parameters. By default, it is set to \code{NULL}.
#' @param s A \code{"numeric"} value indicating the expert's confidence in the prior knowledge. 
#' This argument takes values on the interval \eqn{[0, N]}, where \eqn{N} is the sample size, and is used
#' to synchronize the support of the prior knowledge and the sample.
#' By default, it is \code{NULL}, and must be modified only if prior information is to be 
#' incorporated in the learning process.
#' @param priorData An object of class \code{"data.frame"}, corresponding to the prior information.
#' @param conditionalfunction The output of the internal function \code{learn.tree.Intervals}.
#' @param mm One of the inputs and the output of the recursive internal function \code{"conditional"}.
#' @param scale A \code{"logical"} value indicating whether to standardize the numeric variables to have mean 0 and standard deviation 1.
#' @return The main function \code{conditionalMethod} returns a list with the name of the parents, 
#' the different intervals and the fitted densities
#' @details The main function, \code{conditionalMethod()}, fits truncated basis functions for the conditioned variable 
#' for each configuration of splits of the parent variables. The domain of the parent variables is splitted
#' in different intervals and univariate functions are fitted in these 
#' ranges. The remaining above described functions are internal to the main function.
#' @seealso \link{printConditional}
#' @examples
#' ## Dataset
#' X <- rnorm(1000)
#' Y <- rbeta(1000, shape1 = abs(X)/2, shape2 = abs(X)/2)
#' Z <- rnorm(1000, mean = Y)
#' data <- data.frame(X = X, Y = Y, Z = Z)
#' 
#' ## Conditional Method
#' parents <- c("X","Y")
#' child <- "Z"
#' intervals <- 2
#'
#' potential <- "MTE"
#' fMTE <- conditionalMethod(data, nameParents = parents, nameChild = child, 
#' numIntervals = intervals, POTENTIAL_TYPE = potential)
#' printConditional(fMTE)
#' 
#' ##############################################################################
#' \donttest{
#' potential <- "MOP"
#' fMOP <- conditionalMethod(data, nameParents = parents, nameChild = child,
#' numIntervals = intervals, POTENTIAL_TYPE = potential, maxParam = 15)
#' printConditional(fMOP)
#' }
#' ##############################################################################
#' 
#' ##############################################################################
#' ## Internal functions: Not needed to run #####################################
#' ##############################################################################
#' \donttest{
#' domainP <- lapply(parents, function(i)range(data[,i]))
#' names(domainP) = parents
#' domainC <- range(data[, child])
#' t <- conditional(data, nameParents = parents, nameChild = child,
#' domainParents = domainP, domainChild = domainC, numIntervals = intervals,
#' mm = NULL, POTENTIAL_TYPE = potential)
#' printConditional(t)
#' selection <- select(data, nameParents = parents, nameChild = child,
#' domainParents = domainP, domainChild = domainC, numIntervals = intervals,
#' POTENTIAL_TYPE = potential)
#' parent1 <- selection$parent; parent1
#' domainParent1 <- range(data[,parent1])
#' treeParent1 <- learn.tree.Intervals(data, nameParents = parent1,
#' nameChild = child, domainParents = domainParent1, domainChild = domainC,
#' numIntervals = intervals, POTENTIAL_TYPE = potential)
#' BICscoreMoTBF(treeParent1, data, nameParents = parent1, nameChild = child)
#' }
#' 
#' ###############################################################################
#' ###############################################################################
#' 
#' @export
conditionalMethod <- function(data, nameParents, nameChild, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL, scale = FALSE)
{
  if((POTENTIAL_TYPE=="MOP")||(POTENTIAL_TYPE=="MTE")){
    data <- newData(data, nameChild, nameParents)
    # browser()
    ## Domains
    if(is.factor(data[,nameChild])) domainChild <- levels(data[,nameChild])
    else domainChild <- range(data[,nameChild])
    
    domainParents <- lapply(nameParents, 
                            function(pa){
                              if(is.numeric(data[,pa])){
                                return(range(data[,pa]))
                              }else{
                                return(levels(data[,pa])) 
                              } 
                            })
    names(domainParents) <- nameParents
    
    ## Recursive process
    mm <- c()
    mm <- conditional(data, nameParents, nameChild, domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData, scale = scale)
    return(mm)
  }else{
    return(message("Unknown method, please use MOP or MTE"))
  }
}

#'@rdname conditionalmotbf.learning
#'@export
conditional <- function(data, nameParents, nameChild, domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL, scale = FALSE){  
  
  ## select the parent who get the best BIC score when its domain is splitted
  # browser()
  f <- select(data, nameParents, nameChild, domainChild, domainParents, numIntervals, POTENTIAL_TYPE, maxParam, s, priorData, scale = scale)
  nameParents <- nameParents[which(nameParents!=f$parent)]
  
  for(i in 1:length(f$t)){
    m <- list(parent=f$parent, interval=f$t[[i]]$interval, Px=f$t[[i]]$Px)
    mm[[length(mm)+1]] <- m
    if(is.numeric(data[,f$parent])){
      dataInterval <- splitdata(data, f$parent,
                                f$t[[i]]$interval[1]-0.001*(f$t[[i]]$interval[1]==min(domainParents[[f$parent]])),
                                f$t[[i]]$interval[2])
      if(nrow(dataInterval)>=0||length(nameParents)==0){
        if(length(nameParents)==0){
          mm <- mm
        } else {
          mm[[length(mm)]]$Px <- NULL
          
          ## An recursive process
          mm <- conditional(dataInterval, nameParents, nameChild,domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData, scale = scale)
        }
      }
    } else{
      dataInterval <- subset(data,(data[,f$parent]==f$t[[i]]$interval))
      
      if(nrow(dataInterval)>=0||(length(nameParents)==0)){
        if(length(nameParents)==0){
          mm <- mm
        } else {
          mm[[length(mm)]]$Px <- NULL
          mm <- conditional(dataInterval, nameParents, nameChild,domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData, scale = scale)
        }
      }
    }
  }
  return(mm)
}


#'@rdname conditionalmotbf.learning
#'@export
select <- function(data, nameParents, nameChild, domainChild, domainParents, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL, scale=FALSE)
{
  # browser()
  bestbic <- -10^10; bestvalues <- 0; Bic <- 0
  for(i in 1:length(nameParents)){
    t <- learn.tree.Intervals(data, nameParents[i], nameChild, domainParents, domainChild, numIntervals, POTENTIAL_TYPE, maxParam, s, priorData, scale = scale)
    
    Bic <- BICscoreMoTBF(t, data, nameParents[i], nameChild)
    values <- list(parent=nameParents[i], t=t)
    
    if(bestbic<Bic){
      bestbic <- Bic
      bestvalues <- values
    }
  }
  return(bestvalues)      
}

#'@rdname conditionalmotbf.learning
#'@export
learn.tree.Intervals <- function(data, nameParents, nameChild, domainParents, domainChild, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL, scale=FALSE)
{
  
  X <- data[, nameParents]
  Y <- data[, nameChild]
  if(is.numeric(domainParents)) Xrange <- domainParents
  else Xrange <- domainParents[nameParents][[1]]
  
  NY <- length(Y)
  if(is.numeric(X)){
    priorD <- priorData[,nameChild]
    bestb <- min(X)-1
    B <- quantileIntervals(X, numIntervals)
    Pp1 <- list(); points <- c()
    for(i in 1:length(B)){
      Yf <- Y[which(X>bestb)]
      Xf <- X[which(X>bestb)]
      if(length(Y)==0) next
      if(is.factor(Y)){
        
        # cutPoints <- discreteVariablesStates(nameChild, data)[[1]]$states
        cutPoints = levels(Y)
        # pD <- probDiscreteVariable(cutPoints, Yf)
        pD <- probDiscreteVariable(Yf)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        prob <- list(pD)
        bestBIC <-  getBICDiscreteBN(prob)
      } else{
        
        if(is.null(priorData)) P <- univMoTBF(Yf, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
        # else P <- learnMoTBFpriorInformation(priorD, Yf, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        else P <- learnMoTBFpriorInformation(priorD, Yf, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
        bestBIC <- BICMoTBF(P, Yf)
      }
      b <- B[i]
      Xl <- Xf[Xf<=b]; Xl <- sort(Xl)
      Yl <- Yf[which(Xf<=b)]
      if(length(Yl)<=5) break
      if(length(Yl)==0) next
      # if(max(Yl)==min(Yl)) next # max() and min() are not meaningful for factors
      
      ## discrete child variable
      if(is.factor(Yl)){
        if(length(unique(Yl))==1) next
        cutPoints <- unique(Yl)
        # pD <- probDiscreteVariable(cutPoints, Yl)
        pD <- probDiscreteVariable(Yl)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        Px1 <- pD
      } else{
        
        if(max(Yl)==min(Yl)) next
        
        if(is.null(priorD)) Px1 <- univMoTBF(Yl, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale) 
        # else Px1 <- learnMoTBFpriorInformation(priorD, Yl, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        else Px1 <- learnMoTBFpriorInformation(priorD, Yl, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
      }
      
      Xr <- Xf[Xf>b]; Xr <- sort(Xr)
      Y22 <- Yf[which(Xf>b)]; Y22 <- sort(Y22)
      if(length(Y22)==0) next
      # if(max(Y22)==min(Y22)) next
      
      if(is.factor(Y22)){
        
        if(length(unique(Y22))==1) next
        
        cutPoints <- unique(Y22)
        # pD <- probDiscreteVariable(cutPoints, Y22)
        pD <- probDiscreteVariable(Y22)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        Px2 <- pD
      } else{
        
        if(max(Y22)==min(Y22)) next
        
        if(is.null(priorD)) Px2 <- univMoTBF(Y22, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
        # else Px2 <- learnMoTBFpriorInformation(priorD, Y22, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        else Px2 <- learnMoTBFpriorInformation(priorD, Y22, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
      }
      
      if(is.factor(Y)){
        DiscreteBN <-  list(Px1, Px2)
        BICT <-  getBICDiscreteBN(DiscreteBN)
      } else {
        PX <- list(Px1,Px2)
        multiY <- list(Yl, Y22)
        BICT <- BICMultiFunctions(PX, multiY)
      }
      
      if(is.na(BICT)) next ###quitar
      
      if(BICT>bestBIC){
        Pp1[[length(Pp1)+1]] <- Px1
        Pp2 <- Px2
        bestb <- b
        points <- c(points,bestb)
      }
    }
    
    if(is.null(points)){
      Y <- Y[which(X==X)]
      if(is.factor(Y)){
        cutPoints <- levels(Y)
        # pD <- probDiscreteVariable(cutPoints, Y)
        pD <- probDiscreteVariable(Y)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        prob <- list(pD)
        v <- c()
        values <- list(Px=prob, interval=c(Xrange[1], Xrange[2]))
      } else{
        v <- c()
        if(is.null(priorData)) P <- univMoTBF(Y, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
        # else P <- learnMoTBFpriorInformation(priorD, Y, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        else P <- learnMoTBFpriorInformation(priorD, Y, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
        values <- list(Px=P, interval=c(Xrange[1], Xrange[2]))
      }
      v[[length(v)+1]] <- values
      return(v)
    } else{
      if(length(points)==1){
        v <- c(); X <- data[,nameParents]; X <- X[X<=points]
        values <- list(Px=Pp1[[1]], interval=c(Xrange[1], points))
        v[[length(v)+1]] <- values

        X <- data[, nameParents]; X <- X[X>points]
        values <- list(Px=Pp2, interval=c(points, Xrange[2]))
        v[[length(v)+1]] <- values
      } else{
        v <- c(); X <- data[,nameParents]
        dom <- c(Xrange[1], points, Xrange[2])
        Pp <- Pp1; Pp[[length(Pp)+1]] <- Pp2
        for(i in 1:length(Pp)){
          values <- list(Px=Pp[[i]], interval=c(dom[i], dom[i+1]))
          v[[length(v)+1]] <- values
        }
      }
      
      return(v) 
    }
  } else {# Discrete parent
    
    # B <- discreteVariablesStates(nameParents, data)[[1]]$states
    B = levels(data[,nameParents])
    v <- c()
    for(i in 1:length(B)){
      # browser()
      Y1 <- Y[which(X==B[i])]
      if(is.factor(Y1)){
        cutPoints <- levels(Y1)
        # pD <- probDiscreteVariable(cutPoints, Y1)
        pD <- probDiscreteVariable(Y1)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        values <- list(Px=pD, interval=B[i])
        v[[length(v)+1]] <- values
        
      } else{
        if(length(Y1)==0){# When a combination of parents is not observed
          # assign a uniform distribution
          # browser()
          # P <- asMOPString(1/diff(domainChild))
          P <- do.call(paste0("as",POTENTIAL_TYPE,"String"), list(1/diff(domainChild)))
          P <- list(Function = P, Subclass = tolower(POTENTIAL_TYPE), Domain = domainChild)
          # P <- motbf(P)
          # class(P) <-  c(class(P), tolower(POTENTIAL_TYPE))
          P = do.call(paste0("new_",tolower(POTENTIAL_TYPE)), list(P))
          # next
        }else{
          if(is.null(priorData)){
            P <- univMoTBF(Y1, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale) 
          }else{
            priorChild <- priorData[,nameChild]
            if(ncol(priorData)!=1){
              priorParent <- priorData[,nameParents]
              priorD <- priorChild[which(priorParent==B[i])]
              # P <- learnMoTBFpriorInformation(priorD, Y1, s, POTENTIAL_TYPE, domainChild)$posteriorFunction
              P <- learnMoTBFpriorInformation(priorD, Y1, s, POTENTIAL_TYPE, domainChild, scale = scale)
            } else {
              P <- univMoTBF(Y1, POTENTIAL_TYPE, domainChild, maxParam=maxParam, scale = scale)
            }
          }
        }

        
        values <- list(Px=P, interval=B[i])
        v[[length(v)+1]] <- values
      } 
    }
    return(v)
  }
}

#'@rdname conditionalmotbf.learning
#'@export
BICscoreMoTBF <- function(conditionalfunction, data, nameParents, nameChild)
{
  X <- data[, nameParents]
  Y <- data[, nameChild]
  
  ## Compute the density values
  valuesT <- c(); DiscreteBN <- c(); nlevel <- c()
  for(i in 1:length(conditionalfunction)){
    if(is.numeric(X)){
      
      ##Continuous parent
      if(is.numeric(Y)){
        
        ## Continuous child
        domain <- Y[which((X>=conditionalfunction[[i]]$interval[1])&(X<=conditionalfunction[[i]]$interval[2]))]
        if(is.motbf(conditionalfunction[[i]]$Px)) values <- as.function(conditionalfunction[[i]]$Px)(domain)
        else values <- NULL
        if((length(values)!=length(domain))&&(is.numeric(conditionalfunction[[i]]$Px))) values <- rep(values, length(domain))
        valuesT <- c(valuesT,values)
      } else{
        
        ## Discrete child
        if(length(conditionalfunction[[i]]$Px)==1) DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px[[i]]
        else DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px
      } 
    } else 
      
      ## Discrete Parent
      if(is.numeric(Y)){ 
        # browser()
        ## Continuous child
        domain <- Y[which(X==conditionalfunction[[i]]$interval)] 
        if(is.motbf(conditionalfunction[[i]]$Px)) values <- as.function(conditionalfunction[[i]]$Px)(domain)
        else values <- NULL
        if((length(values)!=length(domain))&&(is.numeric(conditionalfunction[[i]]$Px))) values <- rep(values, length(domain))
        valuesT <- c(valuesT,values)
      } else {
        
        ## Discrete child
        if(length(conditionalfunction[[i]]$Px)==1) DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px[[i]]
        else DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px
      }
  }
  
  ## Compute the BIC
  if(is.numeric(X)){
    if(is.numeric(Y)){
      sizePx=c()
      for(i in 1:length(conditionalfunction)){
        if(!is.motbf(conditionalfunction[[i]]$Px)) tam <- 0
        else tam <- length(coef(conditionalfunction[[i]]$Px))
        sizePx <- c(sizePx, tam)
      }
      sizeP <- sum(sizePx)
      BiC <- sum(log(valuesT))-1/2*(sizeP+(length(conditionalfunction)-1))*log(length(valuesT))
    } else {
      BiC <- getBICDiscreteBN(DiscreteBN)
    }
  } else {
    if(is.numeric(Y)){
      sizePx=c()
      for(i in 1:length(conditionalfunction)){
        if(!is.motbf(conditionalfunction[[i]]$Px)) tam <- 0
        else tam <- length(coef(conditionalfunction[[i]]$Px))
        sizePx <- c(sizePx, tam)
      }
      sizeP <- sum(sizePx)
      valuesT <- valuesT[valuesT>0]
      BiC <- sum(log(valuesT))-1/2*(sizeP+(length(conditionalfunction)-1))*log(length(valuesT))
    } else {
      BiC <- getBICDiscreteBN(DiscreteBN)
    }
  }
  return(BiC)
}


#' BIC score for multiple functions
#' 
#' Compute the BIC score using more than one probability functions.
#' 
#' @param Px A list of objects of class \code{"motbf"}.
#' @param X A list with as many \code{"numeric"} vectors as densities in \code{Px},
#' used to compute the BIC score for each density.
#' @return The \code{"numeric"} BIC value.
#' @seealso \link{univMoTBF}
#' @export
#' @examples
#' ## Data
#' X <- rnorm(500)
#' Y <- rnorm(500, mean=1)
#' data <- data.frame(X=X, Y=Y)
#' ## Data as a "list"
#' Xlist <- sapply(data, list)
#' 
#' ## Learning as a "list"
#' Plist <- lapply(data, univMoTBF, POTENTIAL_TYPE="MOP")
#' Plist
#' 
#' ## BIC value
#' BICMultiFunctions(Px=Plist, X=Xlist)
#' 
BICMultiFunctions <- function(Px, X){
  coeffs <- unlist(sapply(1:length(Px), function(i) coef(Px[[i]])))
  totalValues <- unlist(sapply(1:length(Px), function(i) as.function(Px[[i]])(X[[i]])))
  BiC <- sum(log(totalValues))-(1/2*((length(coeffs)+length(Px))*log(length(totalValues))))
  return(BiC)
}


#' Summary of conditional MoTBF densities
#' 
#' Print the description of an MoTBF demnsity for one variable conditional on another variable.
#' 
#' @param conditionalFunction the output of the function \code{conditionalMethod}. A list with the interval of the parent
#' and the final \code{"motbf"} density function fitted in each interval.
#' @return The results  of the conditional function are shown.
#' @seealso \link{conditionalMethod}
#' @export
#' @examples
#' ## Data
#' X <- rexp(500)
#' Y <- rnorm(500, mean=X)
#' data <- data.frame(X=X,Y=Y)
#' cov(data)
#' 
#' ## Conditional Learning
#' parent <- "X"
#' child <- "Y"
#' intervals <- 5
#' potential <- "MOP"
#' P <- conditionalMethod(data, nameParents=parent, nameChild=child, 
#' numIntervals=intervals, POTENTIAL_TYPE=potential)
#' printConditional(P)
#' 
printConditional <- function(conditionalFunction)
{
  for(i in 1:length(conditionalFunction)){
    if((conditionalFunction[[i]]$parent)==(conditionalFunction[[1]]$parent)){
      if(length(conditionalFunction[[i]])==3){
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
        print(              conditionalFunction[[i]]$Px)
      } else{
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
      }
    }else{
      if(length(conditionalFunction[[i]])==3){
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
        print(              conditionalFunction[[i]]$Px)
      } else{
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
      }
    }
  }
}
