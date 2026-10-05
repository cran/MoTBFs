#' Learning hybrid BNs with MoTBFs
#'
#' Learn mixtures of truncated basis functions in a full hybrid network.
#' 
#' @param graph A network of the class \code{"bn"} (from bnlearn package).
#' @param data An object of class \code{"data.frame"}; it can contain continuous and discrete variables.
#' @param numIntervals A positive integer indicating the maximum number of intervals 
#' for splitting the domain of the continuous parent variables. By default, it is set to 4.
#' @param POTENTIAL_TYPE A \code{"character"} string, either \emph{MOP} or \emph{MTE}, 
#' corresponding to the type of basis function. By default, it is set to \emph{MOP}.
#' @param maxParam A positive integer which indicates the maximum number of coefficients in the function. 
#' If specified, the output is the function which gets the best BIC with, at most, this number of parameters.
#' By default, it is set to \code{NULL}.
#' @param s A \code{"numeric"} value which specifies the expert confidence in the prior knowledge. 
#' This argument takes values on the interval \eqn{[0, N]}, where \eqn{N} is the sample size, and is used
#' to synchronize the support of the prior knowledge and the sample.
#' By default, it is \code{NULL}, and must be modified only if prior information is to be incorporated to the fits.
#' @param priorData An object of class \code{"data.frame"}, corresponding to the prior information.
#' @param scale A \code{"logical"} value indicating whether to standardize the numeric variables to have mean 0 and standard deviation 1.
#' @return A list of lists. Each list contains two elements
#' \item{Child}{A \code{"charater"} string which contains the name of the child variable.}
#' \item{functions}{A list of three elements: the name of the parents; a \code{"numeric"} vector
#' indicating the interval of the parent; and the fitted function in this interval.}
#' @details If the variable is discrete then it computes the probabilities and the size of each leaf. 
#' Children that have discrete parents have as many functions as configurations of the parents. 
#' Children that have continuous parents have as many functions as the number indicated in the 
#' argument \code{"numIntervals"} for each parent. Children that have mixed parents, combine both methods.
#' The BIC criterion is used to decide the number of splitting points of the parent domains and to choose
#' the number of basis functions used.
#' @export
#' @examples
#' 
#' ## Dataset Ecoli
#' require(MoTBFs)
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' intervals <- 3
#' potential <- "MOP"
#' bn1 <- motbf.fit(graph = dag, data = data, numIntervals = intervals, 
#' POTENTIAL_TYPE = potential, maxParam = 5)
#' bn1
#' 
#' ## Learning BN
#' intervals <- 4
#' potential <- "MTE"
#' bn2 <- motbf.fit(graph = dag, data = data, numIntervals = intervals, 
#' POTENTIAL_TYPE = potential, maxParam = 15)
#' bn2
#' 
#' 
motbf.fit <- function(graph, data, numIntervals = 4, POTENTIAL_TYPE = 'MOP', 
                            maxParam = NULL, s = NULL, priorData = NULL, scale = TRUE){

  # browser()
  data = check_data(data)
  if(!is.null(priorData)){
    priorData <- check_data(priorData)
    if(scale){
      priorData_original <- priorData
      priorData <- rescale_data(priorData)
    }
  }
  if(scale){
    data_original = data
    data = rescale_data(data)
  }

  childrenAndParents <- getChildParentsFromGraph(graph,colnames(data))
  MoTBFs <- c()
  for(i in 1:length(childrenAndParents)){

    message(" Learning variable",childrenAndParents[[i]][1], '\n')
    
    if(length(childrenAndParents[[i]])==1){
      
      ## Single variable
      Child <- childrenAndParents[[i]]
      if(Child%in%colnames(priorData)) priorParent <- priorData[,Child]
      else priorParent <- NULL
      
      if(is.numeric(data[,Child])){
        
        ## Numeric child
        if(is.null(priorParent)) PX <- univMoTBF(data[,Child], POTENTIAL_TYPE, maxParam=maxParam, scale = FALSE)
        # else PX <- learnMoTBFpriorInformation(priorParent, data[,Child], s, POTENTIAL_TYPE, maxParam=maxParam)$posteriorFunction 
        else PX <- learnMoTBFpriorInformation(priorParent, data[,Child], s, POTENTIAL_TYPE, maxParam=maxParam, scale = FALSE)
        UnivMoTBFs <- list(PX)
        information <- list(Child=Child, functions=UnivMoTBFs, varType = "Continuous")
        MoTBFs[[length(MoTBFs)+1]] <- information
      }else{
        
        ##Discrete child
        # states <- discreteVariablesStates(Child, data)[[1]]$states
        # prob <- probDiscreteVariable(states, data[,Child])
        prob <- probDiscreteVariable(data[,Child])
        UnivDisc <- list(prob)
        information <- list(Child=Child, functions=UnivDisc, varType = "Discrete")
        MoTBFs[[length(MoTBFs)+1]] <- information
      }
    } else{
      
      ## Conditional variables
      Child <- childrenAndParents[[i]][1]
      Parents <- childrenAndParents[[i]][2:length(childrenAndParents[[i]])]
      if((Child%in%colnames(priorData))&&(all(Parents%in%colnames(priorData)))){
        priorChild <- priorData[,childrenAndParents[[i]]]
        colnames(priorChild) <- childrenAndParents[[i]]
      } else if (Child%in%colnames(priorData)) {
        priorChild <- as.data.frame(priorData[,Child])
        colnames(priorChild) <- Child
      } else{ 
        priorChild <- NULL
      }
      varType <- ifelse(is.numeric(data[,Child]), "Continuous", "Discrete")
      result <- conditionalMethod(data, Parents, Child, numIntervals, POTENTIAL_TYPE, maxParam=maxParam, s, priorChild, scale = FALSE)
      information <- list(Child=Child, functions=result, varType = varType)
      MoTBFs[[length(MoTBFs)+1]] <- information
    }
  }
  
  bn = getFormatedBN(MoTBFs)
  #browser()

  if(scale){
    bn = rescale_motbf_fit(bn, POTENTIAL_TYPE = POTENTIAL_TYPE, data = data_original)
  }
  return(bn)
}

#' @rdname motbf.fit
#' @export
MoTBFs_Learning <- motbf.fit

#' Retrieve DAG from BN
#' 
#' Get the Directed acyclic graph (DAG) from a Bayesian network (BN)
#' @param bn An object of class \code{motbf_fit} or \code{bn.fit} (from bnlearn package).
#' @return An object of class \code{bn} (same class as in bnlearn package).
#' @importFrom bnlearn model2network
#' @export
#' 
getDAG = function(bn){
  nm = names(bn)
  m = matrix(0,ncol = length(nm), nrow = length(nm), dimnames = list(nm, nm))

  for(i in 1:length(bn)){
    node = bn[[nm[i]]]
    child = node$node
    parent = node$parents
    
    r = which(rownames(m) %in% parent)
    c = which(colnames(m) %in% child)
    
    m[r,c] = 1
  }
  
  st = c()
  for(i in 1:ncol(m)){
    node = colnames(m)[i]
    parents = colnames(m)[which(m[,node] ==1)]
    
    if(length(parents)==0){
      parents = ''
    }else{
      parents = paste0('|',paste0(parents, collapse = ':'))
    }
    
    st = paste0(st,paste0('[',node,parents,']'), collapse = '')
  }
  
  dag = model2network(st, ordering = nm)
  return(dag)
}


#' BIC of a hybrid BN
#' 
#' Compute the BIC score and the loglikelihood from the fitted MoTBFs functions 
#' in a hybrid Bayesian network, i.e., from objects of class motbf.fit.
#' 
#' @name goodnessMoTBFBN
#' @rdname goodnessMoTBFBN
#' @param object The output of the 'motbf.fit()' function
#' @param data The dataset of class \code{data.frame}.
#' @return A numeric value giving the log-likelihood of the BN.
#' @seealso \link{MoTBFs_Learning}
#' @examples
#' 
#' ## Dataset Ecoli
#' require(MoTBFs)
#' data(ecoli)
#' data <- ecoli[,-c(1)] ## remove variable sequence
#' 
#' ## Directed acyclic graph
#' dag <- LearningHC(data)
#' 
#' ## Learning BN
#' intervals <- 3
#' potential <- "MOP"
#' P1 <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' numIntervals = intervals, maxParam = 5)
#' logLikelihood.MoTBFBN(P1, data) ##BIC$LogLikelihood
#' BIC <- BiC.MoTBFBN(P1, data)
#' BIC$BIC
#'
#' ## Learning BN
#' intervals <- 2
#' potential <- "MTE"
#' P2 <- MoTBFs_Learning(graph = dag, data = data, POTENTIAL_TYPE=potential,
#' numIntervals = intervals, maxParam = 10)
#' logLikelihood.MoTBFBN(P2, data) ##BIC$LogLikelihood
#' BIC <- BiC.MoTBFBN(P2, data)
#' BIC$BIC 
#' @export
logLikelihood.MoTBFBN <- function(object, data){
  
  # browser()
  disc <- names(which(sapply(data, is.factor)==TRUE))
  
  if(length(disc)!=0){ 
    states <- discreteVariablesStates(disc, data)
  }
  min_cont = sapply(setdiff(colnames(data),disc),function(vari){
    min(attributes(object)$levels[[vari]])
  },simplify = FALSE)
  
  originalData <- data; valuesT <- c()
  
  for(i in 1:length(object)){
    node = object[[i]]
    fx = node$functions
    n = ncol(fx)
    ## Root nodes
    if(is.null(node$parents)){
      data <- originalData
      Y <- data[, node$node]
      
      if(is.numeric(Y)){
        ## Continuous node
        f = eval(parse(text = paste("f <- function(",node$node,")", fx[[1,n]])))
        values <- f(Y)
        valuesT <- c(valuesT,values)
      } else{
        
        ## Discrete node
        params = fx[,ncol(fx)]
        size = attr(fx[[ncol(fx)]], 'sizeDataLeaf')
        
        st <- sum(size)
        for(j in 1:length(params)){
          values=rep(params[j],size[j]*length(Y)/st)
          valuesT=c(valuesT,values)
        }
      }
    } else{ # nodes with parents
      N <- length(valuesT)
      Parents <- node$parents
      
      fullData <- c(); enter <- c(); data <- originalData
      
      #--------------------------------------------------------------------------#
      
      if(node$type == 'Discrete'){
        # if the node is discrete, the length of the fx data.frame is longer due to its own state values
        # here, we reduce the data.frame so that each row corresponds with a combination of parents and the 
        # probability distribution of the discrete variable is stored as a vector
        
        fx = collapseDiscreteCPD(node)
      }
      
      #--------------------------------------------------------------------------#
      for(j in 1:nrow(fx)){
        data <- originalData
        
        for(k in 1:length(Parents)){
          
          if(is.numeric(originalData[,Parents[k]])){
            lower = fx[[j,Parents[k]]][1]
            if(lower==min_cont[[Parents[k]]]){
              lower = lower-0.001
            }
            upper = fx[[j,Parents[k]]][2]
            data <- splitdata(data, Parents[k] , lower, upper)
          }else{
            state = fx[j, Parents[k]]
            data <- subset(data,(data[,Parents[k]]==state))
          }
        } # end loop Parents
        
        fullData[[length(fullData)+1]] <- data
        enter <- c(enter, "YES")
        
        if(node$type == 'Discrete'){
          coeff = fx[[j, ncol(fx)]]
          sizeDataLeaf = attr(fx[[j, ncol(fx)]], 'sizeDataLeaf')
          
          values <- c(); st <- sum(sizeDataLeaf)
          
          for(m in 1:length(coeff)){
            val <- rep(coeff[m], sizeDataLeaf[m]*length(data[,node$node])/st)
            values <- c(values, val)
          }
        }else{
          Y <- data[,node$node]
          
          if(!is.motbf(node$functions[[j, n]])){
            values <- c()
          }else{
            f = eval(parse(text = paste("f <- function(",node$node,")", node$functions[[j, n]])))
            values <- f(Y)
          }
        }
        
        valuesT <- c(valuesT,values)
        
      }
      
      if((length(valuesT)-N)!=nrow(originalData)) valuesT <- duplicatedValues(fullData, valuesT, N, enter)
    } # end if condition for nodes with parents
  } # end for loop that iterates over each node
  loglike <- sum(log(valuesT))
  return(loglike)
}

#' @rdname goodnessMoTBFBN
#' @export
BiC.MoTBFBN <- function(object, data)
{
  t <- logLikelihood.MoTBFBN(object, data)
  param <- c()
  
  for(i in 1:length(object)){
    node = object[[i]]
    fx = node$functions
    n = ncol(fx)
    for(j in 1:nrow(fx)){
      if(is.motbf(fx[[j, n]])){ 
        param <- c(param,length(coef(fx[[j, n]]))+1)
      }
      else{ 
        param <- c(param,length(fx[[j, n]]))
      }
    }
    
  }
  
  BiC <- t - 1/2*(sum(param))*log(nrow(data)) 
  return(list(LogLikelihood=t, BIC=BiC))
}



collapseDiscreteCPD = function(node){
  fx = node$functions
  n = ncol(fx)
  # if the node is discrete, the length of the fx data.frame is longer due to its own state values
  # here, we reduce the data.frame so that each row corresponds with a combination of parents and the 
  # probability distribution of the discrete variable is stored as a vector
  
  
  # store size of data leaf 
  sizeDataLeaf = attr(fx[[n]], 'sizeDataLeaf')
  
  # find number of states of node
  # statesNode = levels(bn)[[node$node]]
  statesNode = unique(fx[[node$node]])
  nStatesNode = length(statesNode)
  
  # get sequence according to number of states
  s = seq(1, nrow(fx), nStatesNode)
  cpd = fx[,n]
  
  # initialize data.frame and combine with combParents
  CPDs = data.frame(matrix(ncol = 2, dimnames = list(NULL, c(node$node, paste0('CPD_',node$node)))))
  
  if(!is.null(node$parents)){
    # parent node values
    pNodes = fx[1:(n-2)]
    # combParents only changes if the node is discrete. Otherwise, it remains the same as pNodes
    combParents = unique(pNodes)
    CPDs = cbind(combParents, CPDs)
  }
  
  j = 0
  q = 1
  for(q in 1:length(s)){
    j = j+1
    cpd_j = cpd[c(s[q]: (s[q]+nStatesNode-1))]
    attr(cpd_j, 'sizeDataLeaf') = sizeDataLeaf[c(s[q]:( s[q]+nStatesNode-1))]
    CPDs[[j, paste0('CPD_',node$node)]] = list(cpd_j)
    CPDs[[j, paste0('CPD_',node$node)]] = cpd_j
    
    
    CPDs[[j, node$node]] = list(statesNode)
    CPDs[[j, node$node]] = statesNode
  }
  # CPDs contains the probability distribution and the data size as an attribute
  return(CPDs)
  
}

# 
# fun <- function(parents, interval, first, j, pos, states, originalData)
# {
#   pos <- which(parents[(first+1):(j+1)]%in%interval[[first+1]]$parent)
#   
#   if((length(pos)!=1)&&(!is.null(pos))){
#     no <- c()
#     for(t in 1:(length(pos)-1)){
#       if(is.character(interval[first+pos[t]][[1]]$interval)){
#         s <- which(names(which(sapply(data, is.character)==TRUE))==interval[first+pos[t+1]][[1]]$parent)
#         sta <- states[[s]]$states
#         estate1 <- which(sta==interval[first+pos[t]][[1]]$interval)
#         estate2 <- which(sta==interval[first+pos[t+1]][[1]]$interval)
#         if((estate1+1)==estate2) no <- c(no,pos[t], pos[t+1])
#       }else {
#         if(interval[first+pos[t]][[1]]$interval[2]==interval[first+pos[t+1]][[1]]$interval[1]) no <- c(no,pos[t], pos[t+1])
#       }
#     }
#     if(any(duplicated(no))) no <- no[-which(duplicated(no))]
#     if(is.null(no)) f1 <- first+1
#     else f1 <- first+no[length(no)]
#   } else {
#     f1  <- first+1
#   }
#   
#   if(f1==(j+1)) return(first + no[1])
#   if(pos[1]==0) return(0)
#   pos <- fun(parents, interval, f1, j, pos, states, originalData)
# }
# 
# fun1 <- function(pos, interval, first,j, t)
# {
#   if(is.null(interval[[pos-1]]$Px)){
#     return(pos)
#   }else{
#     for(h in (pos-1):first){
#       if(is.null(interval[[h]]$Px)){
#         pos <- h; break
#       } 
#     }
#     if(pos==1) return(pos+1)
#     if(interval[[pos]]$parent!=interval[[j+1]]$parent) return(pos+1)
#     pos <- fun1(pos, interval,first,j, pos)
#   }
# }
# 
# discretePos <- function(pos, interval, sta, estate,first, t)
# {
#   for(h in (pos-1):first){
#     if(interval[[pos]]$parent==interval[[h]]$parent) {
#       if(interval[[h]]$interval!=sta[estate-1]) break
#       else pos <- h
#     }
#   }
#   if(estate==1) return(pos)
#   if(sta[estate-1]==sta[1]) return(pos)
#   if(h==1) return(pos)
#   if(t==pos) {
#     if(is.null(interval[[pos-1]]$Px)){
#       return(pos)
#     } else {
#       sta <- sta[-estate]; pos <- discretePos(pos, interval, sta, estate-1, first, pos)
#     }
#   }
#   
#   pos <- discretePos(pos, interval, sta, estate-1,first, pos)
# }

duplicatedValues <- function(fullData, valuesT, N, enter)
{
  coln=c()
  for(i in 1:length(fullData)){
    if(enter[i]=="NO") next
    coln <- c(coln,as.numeric(rownames(fullData[[i]])))
  }
  ind <- which(duplicated(coln, fromLast=TRUE))
  if((length(ind)!=0)&&(all(ind>=0))) valuesT <- valuesT[-(N+ind)]
  return(valuesT)
}



