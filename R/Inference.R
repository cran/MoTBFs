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
  if(is.character(node)){
   node <-  which(lapply(bn, `[[`, "Child") == node)
  }
  # if(!is.null(bn[[node]]$functions[[1]]$coeff) || !is.null(bn[[node]]$functions[[1]]$Px$coeff))){
  #   discrete <- TRUE
  # }else{
  #   discrete <- FALSE
  # }
  if(bn[[node]]$varType =="Discrete"){
    discrete <-  TRUE
  }else{
    discrete <- FALSE
  }
  return(discrete)
}



#' Get the states of all discrete nodes from a MoTFB-BN
#' 
#' This function returns the states of all discrete node from a list obtained from \link{MoTBFs_Learning}.
#' @param bn A list of lists obtained from \link{MoTBFs_Learning}.
#' @param dag A network of class \code{"bn"}.
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
#'   bn <- MoTBFs_Learning(dag, data, POTENTIAL_TYPE = "MTE")
#' 
#' ## Get the states of the discrete nodes
#' 
#'   discreteStatesFromBN(bn, dag)


discreteStatesFromBN <- function(bn, dag){
  variables <- unique(unlist(lapply(bn, `[[`, "Child")))
  n <- length(variables)
  k <- 0
  
  discreteStates <- list()
  for( i in 1:n){
    if(is.discrete(i, bn)){
      k <- k+1
      if(is.root(variables[i], dag)){
        states <- names(bn[[i]]$functions[[1]]$coeff)
      }else if(!is.null(names(bn[[i]]$functions[[1]]$Px$coeff))){ 
        states <- names(bn[[i]]$functions[[1]]$Px$coeff)
      }else{
        states <- names(bn[[i]]$functions[[1]]$Px[[1]]$coeff)
      }
      
      discreteStates[[k]] <- states
      names(discreteStates)[k] <- variables[i]
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
    root = T
  }else{
    root = F
  }
  return(root)
}


#' Initialize Data Frame
#' 
#' The function \code{r.data.frame()} initializes a data frame with as many columns as nodes in the MoTBF-network. It also asings each column its data type, i.e., numeric or character. In the case of character columns, the states of the variable are extracted from the \code{"bn"} argument and included as levels.
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}.
#' @param dag An object of class \code{"bn"}, representing the graph of the bayesian network.
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
#'   bn <- MoTBFs_Learning(dag, data, POTENTIAL_TYPE = "MTE")
#'   
#' ## Initialize a data.frame containing 3 columns (x, y and z) with their attributes.
#'   r.data.frame(bn, dag)

r.data.frame <- function(bn, dag){
  # Extraer nombre nodos 
  variables <- unique(unlist(lapply(bn, `[[`, "Child")))
  n <- length(variables)
  
  # Crear df vacio
  rdf <- data.frame(matrix(ncol = n, nrow = 0))
  colnames(rdf) <- variables
  
  # Determinar si el nodo es discreto o continuo
  for(i in 1:n){
    if(is.discrete(i, bn) == T){
      
      rdf[,i] <- as.character(rdf[,i])
      
      # AÑADIR ESTADOS COMO ATRIBUTOS
      # encontrar estados variable discreta
      states <- discreteStatesFromBN(bn, dag)
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
    observed = T
  }else{
    observed = F
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

parentValues <- function(node, bn, obs, rdf){
  node_idx <- which(lapply(bn, `[[`, "Child") == node)
  cases <- bn[[node_idx]]$functions
  node_par <- unique(unlist(lapply(cases, `[[`, "parent")))
  
  if(is.null(node_par)){
    return(parent_value = NULL)
  }
  
  parent_sampled_values <- rdf[1,which(colnames(rdf) %in% node_par), drop = F]
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
  
  node_idx <- which(lapply(bn, `[[`, "Child") == node)
  cases <- bn[[node_idx]]$functions
  
  node_par <- unique(unlist(lapply(cases, `[[`, "parent")))
  
  # NODO RAIZ
  if(is.null(node_par)){
    if(is.discrete(node, bn)){
      
      # extraer la probabilidad de cada estado
      fx <- bn[[node_idx]]$functions[[1]]$coeff
      
    }else{
      # extraer funcion de densidad
      fx <- bn[[node_idx]]$functions[[1]]
      
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
    
    # buscar la funcion correspondiente segun valor de los padres
    j = 0
    acotado = NULL
    for(k in 1:length(node_par)){
      int_fx = F
      while (int_fx == F) {
        j = j+1
        if(j > length(cases)){
          acotado = "descartar"
          break
        }
        p <- node_par[k]
        
        parent_value <- evi[1,p]
        
        # casos (listas) en las que p es el padre
        case_p <- which(lapply(cases, `[[`, "parent") == p)
        
        # si 'p' no es el padre del caso 'j', saltar iteracion
        if(cases[[j]]$parent != p){
          next
        }
        
        if(!is.discrete(p, bn)){
          lower_int <- cases[[j]]$interval[1]
          upper_int <- cases[[j]]$interval[2]
          if(lower_int <= parent_value & upper_int >= parent_value){
            int_fx = T
            
            if(is.null(acotado) & length(node_par)>1){
              if(j == max(case_p)){
                acotado <- c(j, length(cases))
              }else{
                acotado <- c(j, case_p[which(case_p == j)+1])
              }
            }
          }
        }else{
          if(cases[[j]]$interval == parent_value){
            int_fx = T
            if(is.null(acotado)  & length(node_par)>1){
              if(j == max(case_p)){
                acotado <- c(j, length(cases))
              }else{
                acotado <- c(j, case_p[which(case_p == j)+1])
              }
            }
          }
        }
      }
    }
    
    # el valor "j" es el indice de la lista donde se encuentra la funcion correcta
    if(!is.null(acotado)){
      if((acotado[1]== "descartar")||(j<acotado[1])|| (j>acotado[2])){
        fx <- NA
        return(fx)
      }
    }
    if(is.discrete(node, bn)){
      if(!is.null(cases[[j]]$Px$coeff)){
        fx <- cases[[j]]$Px$coeff
      }else{
        fx <- cases[[j]]$Px[[1]]$coeff
      }
      
    }
    else{
      fx <- cases[[j]]$Px
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
  n <- length(bn)
  for(i in 1: n){
    if(!is.discrete(bn[[i]]$Child, bn)){
      k <- length(bn[[i]]$functions)
      for(j in 1: k)
        if(!is.null(lapply(bn[[i]]$functions, `[[`, "Px")[[k]])){
          fx <- bn[[i]]$functions[[k]]
          type <- toupper(subclass(fx$Px))
          return(type)
        }
    }
  }
}



#' Generate Samples From Conditional MoTBFs
#' 
#' This function generates a sample from conditional MoTBFs.
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}.
#' @param dag An object of class \code{"bn"}, representing the directed acyclic graph.
#' @param obs A \code{data.frame} containing the observed variables. This argument can be omitted if no variable is observed.
#' @param size A non-negative integer giving the number of instances to be generated.
#' @param force_size \code{logical} indicating if the sample must be of the size indicated. As a default, it is set to TRUE.
#' @return A \code{data.frame} containing the generated sample.
#' @importFrom ggm topOrder
#' @importFrom bnlearn amat
#' @importFrom stats na.omit
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
#' ## Specify the evidence set 
#'   obs <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 1, stringsAsFactors=FALSE)
#'   
#' ## Get the conditional sample
#'   sample_MoTBFs(bn, dag, obs, size = 10)
#'   
sample_MoTBFs<- function(bn, dag, obs = NULL, size, force_size = T){
  
  rdf <- r.data.frame(bn, dag)
  
  
  # obtener orden topologico de las variables en el dag
  topo_idx <- topOrder(amat(dag))
  topo <- colnames(rdf)[topo_idx]
 
  check_size <- 0
  s <- 0
  #for(s in 1:size)
    while(check_size < size){  
      s <- s+1
    # recorrer el vector de nodos ordenado topologicamente
    for(h in 1:length(topo)){
      # identificar nodo
      node <- topo[h]
      node_idx <- topo_idx[h]
      
      # Comprobar que estoy en el nodo correcto
      if(node != bn[[node_idx]]$Child){
        stop("No se ha cogido el nodo correcto")
      }
      
      # Variables muestradas en la iteracion s
      rdf_i <- rdf[s,, drop = F]
      
      # Valor de los padres en la iteracion s
      evi <- parentValues(node, bn, obs, rdf_i)
      
      # Funcion condicionada del nodo
      fx <- findConditional(node, bn, evi)
      
      # Descartar muestra si el valor de los padres es incompatible

      # Muestrear
      if(is.discrete(node, bn)){
        # caso discreto
        if(any(is.na(fx))){
          rdf[s,] <- NA
          
          break
        }
        states <- levels(rdf[,node])
        Y <- sample(states, 1, replace = T, prob = fx)
        
      }else{
        # caso continuo
        if(length(fx)==1){
          rdf[s,] <- NA
          
          break
        }
        Y <- rMoTBF(size = 1, fx = fx)
      }
      
      # guardar resultados
      rdf[s,node] <- Y 
      
    }
    check_size <- nrow(na.omit(rdf))
    if(force_size == F && s == size){
      break
    }
  }  
  rdf[,colnames(obs)] <- obs
  rdf <- na.omit(rdf)
  
  return(rdf)
  
}



#' Forward Sampling
#' 
#' \code{forward_sampling()} returns the conditional distribution of a target variable given a set of oberved variables.
#' The forward sampling algorithm approximates the conditional distribution from a random sample. 
#' @param bn A list of lists obtained from the function \link{MoTBFs_Learning}.
#' @param dag An object of class \code{"bn"}, representing the directed acyclic graph.
#' @param target A character string equal to the name of the variable of interest.
#' @param evi A \code{data.frame} containing the observed variables.
#' @param size A non-negative integer giving the number of instances to be generated.
#' @param force_size \code{logical} indicating if the sample must be of the size indicated. As a default, it is set to TRUE.
#' @param ... Optional arguments passed on to the \code{\link{univMoTBF}} function. \code{evalRange}, \code{nparam} and \code{maxParam} can be specified. \code{POTENTIAL_TYPE} is taken from the 'bn' object.
#' 
#' @references Henrion, M. (1988). Propagating uncertainty in Bayesian networks by probabilistic logic sampling. In Machine Intelligence and Pattern Recognition (Vol. 5, pp. 149-163). North-Holland.
#' @return A list containing the conditional distribution of the target variable and a data.frame with the generated sample.
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
#' ## Specify the evidence set and target variable
#'   obs <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 1, stringsAsFactors=FALSE)
#'   node <- "alm2"
#'   
#' ## Get the conditional distribution of 'node' and the generated sample
#'   forward_sampling(bn, dag, target = node, evi = obs, size = 10, maxParam = 15)
#'   
forward_sampling <- function(bn, dag, target, evi, size, force_size = T,...){
  
  start_time <- Sys.time()
  
  rdf <- sample_MoTBFs(bn, dag, evi, size)
  
  type <- motbf_type(bn)
  
  if(is.discrete(target, bn)){
    states <- unique(rdf[,target])
    var <- rdf[,target]
    fx <- probDiscreteVariable(states, var)
  }else{
    fx <- univMoTBF(rdf[,target], POTENTIAL_TYPE = type, ...)
  }
  
  
  time <- Sys.time() - start_time
  
  message("Processing Time: ", time, "secs\n")
  
  return(list(fx = fx, sample = rdf))
}

