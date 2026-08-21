#' Score-based hybrid Bayesian Network structure learning
#' 
#' Learn the structure of a hybrid Bayesian network using the \bold{hill climbing}
#' local search method. 
#' 
#' @param dataset A dataset with discrete and continuous variables. If the discrete  
#' variables are not of class \code{"factor"}, they are automatically converted. 
#' @param numIntervals A \code{"numeric"} value indicating the number of categories 
#' used when discretizing a continuous variable, corresponding to intervals of
#' equal width. By default it is \code{NULL}, meaning that the continuous variables
#' are not discretized.
#' @details \code{LearningHC()} automatically converts non-numeric variables into factors
#' before calling function \code{hc()} from the \code{bnlearn} package. \code{LearningHC()} can also 
#' be used to discretize the dataset, using the equal width method, before calling \code{hc()}.
#' @return The output is a \code{"bn"} object containing the learned graph.
#' @seealso \link[bnlearn]{hc}
#' @importFrom bnlearn hc
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' ecoli <- ecoli[,-1] ## Sequence Name
#' 
#' ## DAG1
#' dag1 <- LearningHC(ecoli)
#' dag1
#' plot(dag1)
#' 
#' ## DAG2
#' dag2 <- LearningHC(ecoli, numIntervals = 10)
#' dag2
#' plot(dag2)
#' 
#' 
LearningHC <- function(dataset, numIntervals=NULL)
{
  if(is.numeric(numIntervals)){
    ## Get the discrete dataset
    dataset <- discretizeVariablesEWdis(dataset, numIntervals, factor=TRUE)
    
    ## Estimate DAG
    dag <- hc(dataset, score="loglik")
  } else{
    
    ## Get discrete variables as factor
    pos <- which(!sapply(dataset, is.numeric))
    if(length(pos)!=0) for(i in pos) dataset[,i] <- as.factor(dataset[,i])
    
    ## Estimate DAG
    dag <- bnlearn::hc(dataset)
  }
  return(dag)
}

#' Get the list of relations in a graph
#' 
#' Compute the parents of each variable in the graph.
#' 
#' @param graph A directed acyclic graph of the class \code{"graphNEL"},
#' \code{"network"} or \code{"bn"}.
#' @param nameVars A character array containing the names of the variables in the graph. 
#' This parameter is only used when \code{graph} is of class \code{"network"}.
#' @return A list where each element is a vector containing the name of a variable 
#' and its parents in the graph.
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' ecoli <- ecoli[,-1] ## Sequence Name
#' 
#' ## DAG1
#' dag1 <- LearningHC(ecoli)
#' dag1
#' plot(dag1)
#' getChildParentsFromGraph(dag1)
#' 
#' ## DAG2
#' dag2 <- LearningHC(ecoli, numIntervals = 10)
#' dag2
#' plot(dag2)
#' getChildParentsFromGraph(dag2)
#' 
getChildParentsFromGraph <- function(graph, nameVars=NULL){
  graphClass = class(graph)
  
  if("graphNEL" %in% graphClass){
    type <- 1
  }else if("network" %in% graphClass){
    type <- 2
  }else if("bn" %in% graphClass){
    type <- 3
  }else{
    stop('Graph class not supported')
  }
  
  switch(type, 
         
## type <- 1 
{namenodes <- graph@nodes;
 edgenodes <- graph@edgeL; 
 childrenAndParents <- list();
 for(i in 1:length(edgenodes)) childrenAndParents[[length(childrenAndParents)+1]] <- 
   c(namenodes[i],namenodes[edgenodes[[i]]$edges]);
},

## type <- 2
{childrenAndParents = list();
 for(j in 1:length(graph$nodes)){
   Child <- nameVars[graph$nodes[[j]]$idx]
   Parents <- nameVars[graph$nodes[[j]]$parents]
   childrenAndParents[[length(childrenAndParents)+1]] <- c(Child, Parents)
 }
},

## type <- 3
{
  #if(is.null(nameVars)) nameVars <- names(graph$nodes);
 nameVars <- names(graph$nodes);
 childrenAndParents = list();
 for(i in 1:length(graph$nodes)){
   Child <- nameVars[i]
   Parents <- graph$nodes[[i]]$parents
   childrenAndParents[[length(childrenAndParents)+1]] <- c(Child, Parents)
 }
}
 )
return(childrenAndParents)
}

#Función para aprender el dag de un Naive Bayes con la variable target como objetivo
getNB <- function(df, target){
  
  red <- paste("[",target,"]", sep = "")
  children <- colnames(df)[-which(colnames(df) == target)]
  
  for(i in 1:(length(children))){
    nodo <- paste("[",children[i],"|",target, "]", sep= "")
    red <- paste(red, nodo, sep = "")
  }
  
  dag_NB <- model2network(red, ordering = colnames(df))
  return(dag_NB)
}


getTAN = function(data, target){
  dataD = discretizeVariablesEWdis(data, 4, factor = TRUE)
  
  dag = tree.bayes(dataD, target)
  dag = model2network(modelstring(dag), ordering = colnames(data))
  return(dag)
}


#' Hybrid Bayesian Network structure learning
#' 
#' Learn the structure of a hybrid Bayesian network, using a fixed method 
#' (Naive Bayes, NB), a restricted method (Tree augmented Naive Bayes, TAN), 
#' or an unrestricted method (the \bold{hill climbing}, HC, score-based local search method). 
#' 
#' @param data A dataset with discrete and continuous variables. If the discrete  
#' variables are not of class \code{"factor"}, they are automatically converted.
#' @param method A \code{"character"} string indicating the method to learn the structure: 
#' NB (naive Bayes), TAN (Tree augmented Naive Bayes), or HC (hill climbing) are the available options.
#' @param target An optional parameter only used in the case of NB and TAN to specify the class variable.
#' @details \code{getStructure()} automatically converts non-numeric variables into factors
#' before calling function \code{hc()} from the \code{bnlearn} package. In the case of TAN, it converts 
#' all numeric and non-numeric variables into factors (using 4 equal width intervals) 
#' before calling \code{tree.bayes()} from the \code{bnlearn} package.
#' @return The output is a \code{"bn"} object containing the learned graph.
#' @seealso \link[bnlearn]{hc}
#' @importFrom bnlearn hc
#' @importFrom bnlearn tree.bayes
#' @importFrom bnlearn modelstring
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' ecoli <- ecoli[,-1] ## Sequence Name
#' 
#' ## DAG1
#' dag1 <- getStructure(ecoli, method = "HC")
#' dag1
#' plot(dag1)
#' 
#' ## DAG2
#' dag2 <- getStructure(ecoli, method = "TAN", target = "mcg")
#' dag2
#' plot(dag2)
#' 
#' 
getStructure = function(data, method, target = NULL){
  if(method == 'NB'){
    DAG = getNB(data, target)
  }else if(method == 'TAN'){
    DAG = getTAN(data, target)
  }else if(method == 'HC'){
    DAG <- LearningHC(data)
  }else{
    stop("Argument 'method' must be either 'NB', 'TAN' or 'HC'. ")
  }
}