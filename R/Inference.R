# library(bnlearn)
# library(ggm)

#' Check discreteness of a node
#' 
#' This function allows to check whether a node is discrete or not
#' @param node A character (name of node) or numeric (index of node in the bn list) input. 
#' @param bn A list of lists obtained from \link{MoTBFs_Learning}.
#' @return \code{is.discrete} returns TRUE or FALSE depending on whether the node is discrete or not.
#' @export
#' @examples  
#' 
#' ## Create a dataset
#'   # Continuous variables
#'   x <- rnorm(100)
#'   y <- rnorm(100)
#'   
#'   # Discrete variable
#'   z <- sample(letters[1:2],size = 100, replace = TRUE)
#'   
#'   data <- data.frame(C1 = x, C2 = y, D1 = z, stringsAsFactors = FALSE)
#'   
#' ## Get DAG
#'   dag <- LearningHC(data)
#' 
#' ## Learn BN
#'   bn <- MoTBFs_Learning(dag, data, POTENTIAL_TYPE = "MTE")
#' 
#' ## Check wheter a node is discrete or not
#' 
#'   # Using its name
#'   is.discrete("D1", bn)
#'   
#'   # Using its index position
#'   is.discrete(3, bn)  

is.discrete <- function(node, bn) {
  if(bn[[node]]$type =="Discrete"){
    discrete <-  TRUE
  }else{
    discrete <- FALSE
  }
  return(discrete)
}



#' Get the states of all discrete nodes from a MoTFB-BN
#' 
#' This function returns the states of all discrete node from a list obtained from \link{motbf.fit}.
#' @param bn A list of lists obtained from \link{motbf.fit}.
#' @return \code{discreteStatesFromBN} returns a list of length equal to the number of discrete nodes in the network. Each element of the list corresponds to a node and contains a character vector indicating the states of the node.
#' @export
#' @examples 
#' 
#' ## Create a dataset
#'   # Continuous variables
#'   x <- rnorm(100)
#'   y <- rnorm(100)
#'   
#'   # Discrete variable
#'   z <- sample(letters[1:2],size = 100, replace = TRUE)
#'   
#'   data <- data.frame(C1 = x, C2 = y, D1 = z, stringsAsFactors = FALSE)
#'   
#' ## Get DAG
#'   dag <- LearningHC(data)
#' 
#' ## Learn a BN
#'   bn <- motbf.fit(dag, data, POTENTIAL_TYPE = "MTE")
#' 
#' ## Get the states of the discrete nodes
#' 
#'   discreteStatesFromBN(bn)
#'   
discreteStatesFromBN <- function(bn){
  
  variables <- names(bn)
 
  discreteStates <- list()
  k=1
  
  for (i in variables) {
    if(is.discrete(i,bn)){
      discreteStates[[k]]<-unique(bn[[i]]$functions[[i]])
      names(discreteStates)[k]<-i
      k=k+1
    }
  }
  
  return(discreteStates)
}


#' Root nodes
#' 
#' \code{is.root} checks whether a node has parents or not.
#' @param node A character string indicating the node's name.
#' @param dag An object of class \code{"bn"}.
#' @return \code{is.root} returns TRUE or FALSE depending on whether the node is root or not.  
#' @importFrom bnlearn root.nodes
#' @export
#' @examples 
#' 
#' ## Create a dataset
#'   # Continuous variables
#'   x <- rnorm(100)
#'   y <- rnorm(100)
#'   
#'   # Discrete variable
#'   z <- sample(letters[1:2],size = 100, replace = TRUE)
#'   
#'   data <- data.frame(C1 = x, C2 = y, D1 = z, stringsAsFactors = FALSE)
#'   
#' ## Get DAG
#'   dag <- LearningHC(data)
#'   
#' ## Check if a node is root
#'  is.root("C1", dag)

is.root <- function(node, dag){
  r <- root.nodes(dag)
  if(node %in% r){
    root = TRUE
  }else{
    root = FALSE
  }
  return(root)
}


#' Initialize Data Frame
#' 
#' The function \code{r.data.frame()} initializes a data frame with as many columns as nodes in the MoTBF-network. It also asings each column its data type, i.e., numeric or character. In the case of character columns, the states of the variable are extracted from the \code{"bn"} argument and included as levels.
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}.
#' @return An object of class \code{"data.frame"}, which contains the data type of each column and has no rows.
#' @export
#' @examples 
#' 
#' ## Create a dataset
#'   # Continuous variables
#'   x <- rnorm(100)
#'   y <- rnorm(100)
#'   
#'   # Discrete variable
#'   z <- sample(letters[1:2],size = 100, replace = TRUE)
#'   
#'   data <- data.frame(C1 = x, C2 = y, D1 = z, stringsAsFactors = FALSE)
#'   
#' ## Get DAG
#'   dag <- LearningHC(data)
#'   
#' ## Learn a BN
#'   bn <- motbf.fit(dag, data, POTENTIAL_TYPE = "MTE")
#'   
#' ## Initialize a data.frame containing 3 columns (x, y and z) with their attributes.
#'   r.data.frame(bn)

r.data.frame <- function(bn){
  # Extraer nombre nodos 
  variables <- names(bn)
  n <- length(variables)
  
  # Crear df vacio
  rdf <- data.frame(matrix(ncol = n, nrow = 0))
  colnames(rdf) <- variables
  
  # Determinar si el nodo es discreto o continuo
  for(i in 1:n){
    if(is.discrete(i,bn) == TRUE){
      
      rdf[,i] <- as.factor(rdf[,i])
      # AÑADIR ESTADOS COMO ATRIBUTOS
      # encontrar estados variable discreta
      states <- discreteStatesFromBN(bn)
      states_idx <- which(names(states) == colnames(rdf[i]))
      states_node <- states[[states_idx]]
      levels(rdf[,i]) <- states_node
      
    }else{
      rdf[,i] <- as.numeric(rdf[,i])
    }
  }
  return(rdf)
}



#' Observed Node
#' 
#' \code{is.observed()} checks whether a node belongs to the evidence set or not.
#' @param  node A \code{character} string, matching the node's name.
#' @param evi A \code{data.frame} of the evidence set.
#' @return This function returns TRUE if "node" is included in "evi", or, otherwise, FALSE.
#' @export
#' @examples 
#' 
#' ## Data frame of the evidence set
#'   obs <- data.frame(lip = "1", alm2 = 0.5, stringsAsFactors=FALSE)
#'   
#' ## Check if x is in obs
#'   is.observed("x", obs)

is.observed <- function(node, evi){

  if(node %in% colnames(evi)){
    observed = TRUE
  }else{
    observed = FALSE
  }
  return(observed)
}


#' Value of Parent Nodes
#' 
#' This function returns a \code{data.frame} of dimension '1xn' containing the values of the 'n' parents of a 'node' of interest. 
#' Use this function if you have a random sample and an observed sample with information about the parents.
#' The values of the parents are obtained from the evidence set unless they are not observed. In this case, the values are taken from the random sample.
#' @param node A \code{character} string that represents the node's name.
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}. It contains the conditional density functions of the bayesian network.
#' @param obs A \code{data.frame} of dimension '1xm' containing an instance of the 'm' variables that belong to the evidence set.
#' @param rdf A \code{data.frame} of dimension '1xk' containing an instance of the 'k' variables sampled from the bayesian network.
#' @return  A \code{data.frame} containing the values of the parents of 'node'. 
#' @noRd
#' @examples 
#' 
#' ## Dataset
#'   data("ecoli", package = "MoTBFs")
#'   data <- ecoli[,-c(1,9)]
#' 
#' ## Get directed acyclic graph
#'   dag <- LearningHC(data)
#'   
#' ## Learn bayesian network
#'   bn <- MoTBFs_Learning(dag, data = data, numIntervals = 4, POTENTIAL_TYPE = "MTE")
#'   
#' ## Specify the evidence set
#'   obs <- data.frame(lip = "1", alm1 = 0.5, stringsAsFactors=FALSE)
#'   
#' ## Create a random sample
#'   contData <- data[ ,which(lapply(data, is.numeric) == TRUE)]
#'   fx <- lapply(contData, univMoTBF, POTENTIAL_TYPE = "MTE")
#'   disData <- data[ ,which(lapply(data, is.numeric) == FALSE)]
#'   conSample <- lapply(fx, rMoTBF, size = 1)
#'   disSample <- lapply(unique(disData), sample, size = 1)
#'   
#'   rdf <- as.data.frame(list(conSample,disSample), stringsAsFactors = FALSE)
#'   
#' ## Get the values of the parents of node "alm2"
#'   parentValues("alm2", bn, obs, rdf)
#' 
# FUNCION OBSOLETA
parentValues <- function(node, bn, obs, rdf){
 
  node_par <- bn[[node]]$parents
  
  if(is.null(node_par)){
    return(parent_value = NULL)
  }
  
  parent_sampled_values <- rdf[1,which(colnames(rdf) %in% node_par), drop = FALSE]
  parent_value <- parent_sampled_values
  for(i in 1:length(node_par)){
    p <- node_par[i]
    # si el padre esta incluido en el conjunto de variables observadas 'obs', 
    # se sustituye el valor muestreado por el observado
    if(is.observed(p, obs)){
      parent_value[,p] <- obs[1,p]
    }else{
      parent_value[,p] <- parent_sampled_values[1,p]
    }
  }
  
  return(parent_value)
}
  

#' Find Fitted Conditional MoTBFs
#' 
#' This function returns the conditional probability function of a node given an MoTBF-bayesian network and the value of its parents.
#' @param node A \code{character} string, representing the tardet variable.
#' @param bn A list of lists obtained from \link{MoTBFs_Learning}, containing the conditional functions.
#' @param evi A \code{data.frame} of dimension '1xn' that contains the values of the 'n' parents of the target node. 
#' This argument can be \code{NULL} if \code{"node"} is a root node.
#' @return A list containing the conditional distribution of the target variable.
#' @export
#' @examples 
#' 
#' ## Dataset
#'   data("ecoli", package = "MoTBFs")
#'   data <- ecoli[,-c(1,9)]
#' 
#' ## Get directed acyclic graph
#'   dag <- LearningHC(data)
#'   
#' ## Learn bayesian network
#'   bn <- MoTBFs_Learning(dag, data = data, numIntervals = 4, POTENTIAL_TYPE = "MTE")
#'   
#' ## Specify the evidence set and node of interest
#'   evi <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 1, stringsAsFactors=FALSE)
#'   node = "alm2"
#'   
#' ## Get the conditional distribution
#'   findConditional(node, bn, evi)
#' 
findConditional <- function(node, bn, evi = NULL){
  
  if(!is.null(evi)){
    if(nrow(evi)>1){
      stop("There are more than one value in the observation of the parents")
    }
  }
  
  node_par <- bn[[node]]$parents
  
  # NODO RAIZ
  if(is.null(node_par)){
    if(is.discrete(node,bn)){
      
      # extraer la probabilidad de cada estado
      indx<-paste("CPD_",node,sep = "")
      distr<-as.vector(bn[[node]]$functions[[indx]])
      names(distr)<-bn[[node]]$functions[[node]]
      fx <- distr
      
    }else{
      # extraer funcion de densidad
      indx<-paste("CPD_",node,sep = "")
      fx <- bn[[node]]$functions[[indx]][[1]]
    }
    # NODOS HIJOS
  }else{
    ## Check the evidence set
    # if the evidence set is empty
    if(is.null(evi)){
      stop("The evidence set of the parents is not defined")
    }
    
    # comprobar que la lista de padres es correcta
    if(any(!(node_par %in% colnames(evi)))){
      stop(paste("Some parents of node",node,"are not included in the evidence set"))
    }
    
    #For each parent node we check the positions where its observation matches in functions data frame
    coinci<-list()
    cont=1
    h=nrow(bn[[node]]$functions)
    
    for (k in node_par) {
      
      p=evi[1,k]
      adi<-c()
      
      if(is.discrete(k,bn)){
        for(j in 1:h){
          if(p==bn[[node]]$functions[[k]][[j]]){
            adi<-append(adi,j)
          }
        }
        
      }else{
        
        for (j in 1:h) {
          exi=bn[[node]]$functions[[k]][[j]]['min',]
          exd=bn[[node]]$functions[[k]][[j]]['max',]
          
          if(p>=exi & p<=exd){
            adi<-append(adi,j)
          }
        }
      }
      
      coinci[[cont]]=adi
      cont=cont+1
    }
    
    names(coinci)=node_par
    
    #We look for the intersection of all observations and their positions in functions data frame
    
    inter<-coinci[[1]]
    if(length(coinci)==1){
      pos<-min(inter)
    } else{
      for (i in 1:(length(node_par)-1)) {
      inter<-intersect(inter,coinci[[i+1]])
     }
    }
    pos<-min(inter)
    
    #We select the distribution in position pos in functions data frame
    
    indx<-paste("CPD_",node,sep = "")
    
    if(is.discrete(node,bn)){
      sts<-length(discreteStatesFromBN(bn)[[node]])
      distr<-as.vector(bn[[node]]$functions[[indx]][c(pos:(pos+sts-1))])
      names(distr)<-bn[[node]]$functions[[node]][c(1:sts)]
      fx <- distr
    }else{
      fx<-bn[[node]]$functions[[indx]][[pos]]
    }
  }
    return(fx)
}




#' Type of MoTBF
#' 
#' This function checks whether the density functions of a MoTBF-BN are of type MTE or MOP.
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}.
#' @return A character string, specifying the subclass of MoTBF, i.e., either MTE or MOP.
#' @export
#' @examples 
#' 
#' ## Dataset
#'   data("ecoli", package = "MoTBFs")
#'   data <- ecoli[,-c(1,9)]
#' 
#' ## Get directed acyclic graph
#'   dag <- LearningHC(data)
#'   
#' ## Learn bayesian network
#'   bn <- MoTBFs_Learning(dag, data = data, numIntervals = 4, POTENTIAL_TYPE = "MTE") 
#'   
#' ## Get MoTBF sub-class
#'   motbf_type(bn)

motbf_type <- function(bn){
  subclass = unique(sapply(bn, '[[', 'subclass'))
  type = toupper(subclass[subclass!= "Multinomial"])
  if(length(type)==0){
    type = 'MOP'
  }
  # n <- length(bn)
  # for(i in 1: n){
  #   if(!is.discrete(bn[[i]]$node, bn)){
  #     type <- toupper(bn[[i]]$subclass)
  #     break
  #   }
  # }
  return(type)
}





#' Approximate inference
#' 
#' \code{get_approx_posterior()} returns an approximation to the posterior probability distribution 
#' of a target variable given a set of observed variables. The inference process is based on sample generation. See details.
#' 
#' @param bn An object of class \code{motbf_fit}, obtained from function \link{motbf.fit}.
#' @param target A character string equal to the name of the variable of interest.
#' @param evidence A \code{data.frame} of one row containing the value of the observed variables.
#' @param size A non-negative integer giving the number of random samples to generate from \code{bn}.
#' @param parallel \code{logical} indicating if the particle generation should be parallelized. As a default, it is set to FALSE.
#' @param ... Optional arguments passed on to the \code{\link{univMoTBF}} function. \code{evalRange}, \code{nparam} and \code{maxParam} can be specified. \code{POTENTIAL_TYPE} is taken from the 'bn' object.
#' 
#' @details
#' If any node is observed, i.e., argument \code{evidence} is not NULL, 
#' samples are generated from the Bayesian network using the likelihood weighting algorithm. 
#' Otherwise, i.e., no node is observed, samples are generated using the forward sampling algorithm.
#' 
#' @references Henrion, M. (1988). Propagating uncertainty in Bayesian networks by probabilistic logic sampling. In Machine Intelligence and Pattern Recognition (Vol. 5, pp. 149-163). North-Holland.
#' @return A list of two elements: 1) the posterior probability distribution of the target variable, and 
#' 2) a data.frame with the generated sample, whose weights are attached as an attribute called weights (if \code{evidence} is not NULL).
#' @export
#' @examples 
#' 
#' ## Dataset
#'   data("ecoli", package = "MoTBFs")
#'   data <- ecoli[,-c(1,9)]
#' 
#' ## Get directed acyclic graph
#'   dag <- LearningHC(data)
#'   
#' ## Learn bayesian network
#'   bn <- motbf.fit(dag, data = data, numIntervals = 4, POTENTIAL_TYPE = "MOP")
#'   
#' ## Specify the evidence set and target variable
#'   obs <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 1, stringsAsFactors=FALSE)
#'   node <- "alm2"
#'   
#' ## Get the posterior distribution of 'node' given "evidence" and the generated sample
#'   get_approx_posterior(bn, target = node, evidence = obs, size = 10, maxParam = 15)
#'   
get_approx_posterior <- function(bn, target, evidence = NULL, size = 100, parallel = FALSE,...){
  
  start_time <- Sys.time()
  
  evi = evidence
  evi<-evi[which(names(evi)!=target)]
  
  rdf <- sample_motbfs(bn, n = size, evidence = evi, parallel = parallel)
  
  y = rdf[,target]
  if(!is.null(evidence)){
    w = attr(rdf, 'w')
    wn = w/sum(w)
    y = sample(y, size*10, replace = TRUE, prob = wn)
  }

  
  
  type <- motbf_type(bn)
  
 
  
  if(is.discrete(target, bn)){
    states <- levels(bn)[[target]]
    # var <- factor(rdf[,target], levels = states)
    var <- factor(y, levels = states)
    fx <- probDiscreteVariable(var)
  }else{
    fx <- univMoTBF(y, POTENTIAL_TYPE = type, ...)
    # fx <- univMoTBF(rdf[,target], POTENTIAL_TYPE = type, ...)
  }
  
  
  time <- Sys.time() - start_time
  
  message("Processing Time: ", time, "secs\n")
  
  return(list(fx = fx, sample = rdf))
}





one_sample_motbfs = function(bn, sampling_order, rdf){
  
  for(i in 1:length(sampling_order)){
    # identificar nodo
    node <- sampling_order[i]
    
    
    # Variables muestradas en la iteracion s
    rdf_i <- rdf[1,, drop = FALSE]
    
    # Valor de los padres en la iteracion s
    evi = rdf_i[,bn[[node]]$parents, drop = FALSE]
    
    
    # Funcion condicionada del nodo
    fx <-findConditional(node, bn, evi)
    
    # Descartar muestra si el valor de los padres es incompatible
    
    # Muestrear
    if(is.discrete(node, bn)){
      # caso discreto
      if(any(is.na(fx))){
        rdf[1,] <- NA
        
        break
      }
      states <- levels(rdf[,node])
      Y <- sample(states, 1, replace = TRUE, prob = fx)
      
    }else{
      # caso continuo
      if(length(fx)==1| is.null(fx)){
        rdf[1,] <- NA
        
        break
      }
      Y <- rMoTBF(size = 1, fx = fx)
    }
    
    
    # guardar resultados
    rdf[1,node] <- Y 
    
  }
  return(rdf)
}


one_sample_motbfs_evidence = function(bn, sampling_order, rdf, evidence = NULL){
  # browser()
  w = 1
  for(i in 1:length(sampling_order)){
    # identificar nodo
    node <- sampling_order[i]
    
    
    # Variables muestradas en la iteracion s
    rdf_i <- rdf[1,, drop = FALSE]
    
    # Valor de los padres en la iteracion s
    evi = rdf_i[,bn[[node]]$parents, drop = FALSE]
    
    
    # Funcion condicionada del nodo
    fx <-findConditional(node, bn, evi)
    
    # Nodo observado
    if(node %in% colnames(evidence)){
      Y = evidence[node]
      
      if(is.discrete(node, bn)){
        a = unlist(Y)
        w = w * unname(fx[a])
      }else{
        a = unlist(as.function(fx)(Y))
        w = w * unname(a)
      }
    }else{
      # Muestrear
      if(is.discrete(node, bn)){
        # caso discreto
        if(any(is.na(fx))){
          rdf[1,] <- NA
          
          break
        }
        states <- levels(rdf[,node])
        Y <- sample(states, 1, replace = TRUE, prob = fx)
        
      }else{
        # caso continuo
        if(length(fx)==1| is.null(fx)){
          rdf[1,] <- NA
          
          break
        }
        Y <- rMoTBF(size = 1, fx = fx)
      }
    }
    
    # guardar resultados
    rdf[1,node] <- Y 
    
  }
  return(list(sample = rdf, weights = w))
}


#' Generate Samples From an MoTBF Bayesian network
#' 
#' This function generates a sample from an MoTBF Bayesian network.
#' @param bn An object of class motbf_fit, obtained from the function \link{motbf.fit}.
#' @param n A non-negative integer giving the number of instances to be generated.
#' @param parallel A \code{logical} value. If TRUE, parallelization is carried out. As a default, it is set to FALSE
#' @param evidence A data.frame of one row containing the values for the observed variables. As a default, it is NULL.
#' @return A \code{data.frame} containing the generated sample. Is evidence is not NULL, attribute 'w' contains the weight of each sample.
#' @importFrom ggm topOrder
#' @importFrom bnlearn amat
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @export
#' @examples
#'   data("ecoli", package = "MoTBFs")
#'   dat <- ecoli[,-c(1,4,5,9)]
#'   
#'   # Build DAG
#'   dag <- LearningHC(dat)
#'   
#'   # Learn BN parameters
#'   bn = motbf.fit(dag, dat)
#'  
#'  # Get sample from bn
#'  sam = sample_motbfs(bn, 50)
#'  
#'  
sample_motbfs<- function(bn, n, parallel = FALSE, evidence = NULL){
  
  dag = getDAG(bn)
  
  rdf <- r.data.frame(bn)
  
  
  # obtener orden topologico de las variables en el dag
  topo_idx <- topOrder(amat(dag))
  topo <- colnames(rdf)[topo_idx]
  
  if(is.null(evidence)){
    sampling = one_sample_motbfs
    args = list(bn = bn, sampling_order = topo, rdf = rdf)
  }else{
    sampling = one_sample_motbfs_evidence
    args = list(bn = bn, sampling_order = topo, rdf = rdf, evidence = evidence)
  }
  
  if(parallel){
    cores = detectCores()-1
    
    res = mclapply(1:n, function(i){
      do.call(sampling, args)
    }, mc.cores = cores)
    
  }else{
    res = lapply(1:n, function(i){
      do.call(sampling, args)
    })
  }
  
  if(is.null(evidence)){
    rdf <- do.call(rbind, res)
  }else{
    A = lapply(res, '[[',1)
    
    rdf <- do.call(rbind, A)
    w = sapply(res, '[[',2)
    
    attr(rdf, 'weights') = w
  }
  
  
  return(rdf)
  
}
# sample_motbfs<- function(bn, n, parallel = FALSE){
#   
#   dag = getDAG(bn)
#   
#   rdf <- r.data.frame(bn)
#   
#   
#   # obtener orden topologico de las variables en el dag
#   topo_idx <- topOrder(amat(dag))
#   topo <- colnames(rdf)[topo_idx]
#   
#   if(parallel){
#     cores = detectCores()-1
#     
#     res = mclapply(1:n, function(i){
#       one_sample_motbfs(bn, topo, rdf)
#     }, mc.cores = cores)
#     
#   }else{
#     res = lapply(1:n, function(i){one_sample_motbfs(bn, topo, rdf)})
#   }
#   
#   
#   rdf <- do.call(rbind, res)
#   
#   return(rdf)
#   
# }


#' Conditional probability queries
#' 
#' Compute conditional probability queries from a sample.
#' @param sample A \code{data.frame} containing a sample.
#' @param event A expression describing the event of interest
#' @param evidence A expression describing the conditioning evidence. 
#' @return The conditional probability: P(event|evidence).
#' @export
#' @examples
#'   data("ecoli", package = "MoTBFs")
#'   dat <- ecoli[,-c(1,4,5,9)]
#'   
#'   # Build DAG
#'   dag <- LearningHC(dat)
#'   
#'   # Learn BN parameters
#'   bn = motbf.fit(dag, dat)
#'  
#'  # Get sample from bn
#'  sam = sample_motbfs(bn, 50)
#'  
#'  # Compute P(mcg > 0.6 | aac < 0.8 & alm2 < 0.3)
#'  query(sam, event = (mcg >0.6), evidence = (aac<0.8 & alm2 <0.3))
#'  
query = function(sample, event, evidence){
  
  evi = eval(substitute(evidence), sample)
  res = eval(substitute(event & evidence), sample)
  
  cp = sum(res)/sum(evi)
  return(cp)
}
