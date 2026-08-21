get_Formated_BN = function(bn){
  variables = sapply(bn, "[[", "Child")
  
  fBN <- list()
  levels <- list()
  for(j in 1:length(variables)){
    
    node = variables[j]
    node_idx <- which(lapply(bn, `[[`, "Child") == node)
    cases = bn[[node_idx]]$functions
    
    #
    #------ dataframe para guardar intervalo padres -----#
    nombrePadres = unlist(unique(lapply(cases, `[[`, "parent")))
    
    
    padres <- data.frame(matrix(ncol = length(nombrePadres), nrow = 2))
    colnames(padres) = nombrePadres
    #----------------------------------------#
    
    
    #**** extraer funciones (adaptar si es nodo raiz)
    if(length(nombrePadres)== 0){
      b = list(Px = cases[[1]])
      cases = list(b)
      padres = NULL
    }else{
      # las funciones estan en algunas listas
      b =lapply(cases, `[[`, "Px")
    }
    # donde hay funcion
    index_fx = which(sapply(b, is.null)==FALSE)
    #----------------------------------------#
    
    k = 0
    
    functions <- list()
    
    # loop over cases
    for(i in 1:length(b)){
      nP = cases[[i]]$parent
      padres[,nP] = cases[[i]]$interval
      
      
      if(i %in%index_fx){
        fx = cases[[i]]$Px
        if(length(fx)==1){
          fx = unlist(fx, recursive = FALSE)
        }
        # DISCRETE NODE
        if(bn[[j]]$varType == 'Discrete'){
          states = attr(fx$coeff, 'names')
          px = list()
          px$Function = fx$coeff
          px$Domain = matrix(states, dimnames = list(c(paste0("state", 1:length(states))),node))
          px$sizeDataLeaf = fx$sizeDataLeaf
          subclass = 'Multinomial'
          fx = px
          levels[[j]] <- states
          # CONTINUOUS NODE
        }else{
          # fx$Function = gsub("x", node, fx)
          fx$Function = gsub('EXP', 'exp', gsub('x', node ,gsub('exp\\(', 'EXP(', fx)))
          # fx$Domain <- matrix(as.matrix(fx$Domain), dimnames = list(c('min', 'max'),node))
          fx$Domain <- matrix(fx$Domain, dimnames = list(c('min', 'max'),node))
          subclass = fx$Subclass
          levels[[j]] <- fx$Domain 
        }
        res <- list()
        res['parentInterval'] = list(padres)
        res$Fx = fx
        k = k+1
        
        functions[[k]] = res
      }
      
    }
    
    names(functions) = paste0('Interval_',1:length(functions))
    
    fnode <- list()
    fnode$node = node
    fnode['parents'] = list(nombrePadres)
    fnode$type = bn[[node_idx]]$varType
    fnode$subclass = subclass
    fnode$functions  = functions
    attr(fnode, which = "scope") = c(nombrePadres, node)
    
    fBN[[j]] <- fnode
    
    
  }
  
  names(fBN) <- variables
  names(levels) <- variables
  attr(fBN, which = 'levels') <- levels
  return(fBN)
}



# The CPDs are stored in a data.frame
getFormatedBN_df <- function(bn){
  
  for(j in 1:length(bn)){ #loop over each node of the BN
    node = bn[[j]]
    
    
    # Store conditional probability distributions in a data.frame, including parent's domain
    CPDs = data.frame(matrix(ncol = length(node$parents)+2, 
                             dimnames = list(NULL, c(node$parents,node$node, paste0('CPD_',node$node)))))
    
    fx = node$functions
    nodeType = node$type
    domain = levels(bn)[[node$node]]
    cpt = c()
    sizeDataLeaf = c()
    for(i in 1:length(fx)){# depends on the parent splits
      parentIntervals = fx[[i]]$parentInterval
      CPDs
      
      # typeParent = sapply(parentIntervals, class)
      
      # Fill in the parent's domain
      if(!is.null(parentIntervals)){
        id = match(colnames(CPDs), colnames(parentIntervals))
        id = id[-which(is.na(id))] # remove last column (target node name)
        # browser()
        parentIntervals = parentIntervals[,id, drop = FALSE] # sort "parentIntervals" as "CPDs"
        
        for(k in 1:ncol(parentIntervals)){
          if(is.numeric(parentIntervals[,k])){
            A = parentIntervals[,k, drop = TRUE]
            A = matrix(A, dimnames = list(c('min', 'max'),node$parents[k]))
            CPDs[[i,k]] = list(A)
            CPDs[[i,k]]  = A
          }else{
            CPDs[i, k] = as.character(unique(parentIntervals[,k]))
          }
        }
      }
      
      # target
      if(nodeType=='Discrete'){
        cpt = c(cpt, fx[[i]]$Fx$Function)
        sizeDataLeaf = c(sizeDataLeaf, fx[[i]]$Fx$sizeDataLeaf)
        attr(cpt, 'sizeDataLeaf') = sizeDataLeaf
        
      }else{
        # attributes(domain) = NULL
        CPDs[[i, node$node]] = list(domain)
        CPDs[[i, node$node]] = domain
        
        B = fx[[i]]$Fx
        CPDs[[i, paste0('CPD_',node$node)]] = list(B)
        CPDs[[i, paste0('CPD_',node$node)]] = B
      }
    }
    # target
    
    if(nodeType=='Discrete'){
      CPDs2 = data.frame(sapply(CPDs, rep, each = length(domain))) # converts character vector of discrete parents to list()
      
      disc.parents = which(sapply(CPDs, is.character)) # reconvert to vector
      if(length(disc.parents>0)){
        for(z in 1:length(disc.parents)){
          CPDs2[,disc.parents[z]] = unlist(CPDs2[,disc.parents[z]])
        }
      }
      
      
      CPDs2[,node$node] = names(cpt)
      CPDs2[, paste0('CPD_',node$node)] = cpt
      CPDs = CPDs2
    }
    
    
    bn[[j]]$functions = CPDs
  }
  
  # class(bn) = 'motbf.fit'
  bn = new_motbf_fit(bn)
  return(bn)
}


getFormatedBN <- function(bn){
  fbn = get_Formated_BN(bn)
  dag = getDAG(fbn)
  
  fBN <- list()
  for(i in 1:length(fbn)){
    node = fbn[[i]]
    children = bnlearn::children(dag, node$node)
    if(length(children)==0){
      children = NULL
    }
    
    fnode <- list()
    fnode$node = node$node
    fnode['parents'] = list(node$parents)
    fnode['children'] = list(children)
    fnode$type = node$type
    fnode$subclass = node$subclass
    fnode$functions  = node$functions
    attr(fnode, which = "scope") = c(node$parents, node$node)
    attr(fnode, which = "levels") = attr(fbn, which = 'levels')[c(node$parents, node$node)]
    class(fnode) = "motbf.fit.node"
    fBN[[i]] <- fnode
  }
  names(fBN) <- names(fbn)
  attr(fBN, which = 'levels') <- attr(fbn, which = 'levels')
  
  fBN = getFormatedBN_df(fBN)
  return(fBN)
}