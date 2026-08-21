# filter dataset ('data') to return the cases that match the value or range ('range') of the specified variables ('varName')
filterdata <- function(data, varName, range){
  
  for(i in 1:length(varName)){
    parent <- data[, varName[i], drop = TRUE]
    if(is.numeric(parent)){
      minGlobal = min(parent)
    }
    
    
    dom = unlist(range[i])
    
    if(length(dom)==1){
      data = subset(data,(parent == dom))
      
    }else{
      min = dom[1];max = dom[2]
      
      if(minGlobal == min){
        min = min-1
      }
      data = subset(data,((parent>min)&(parent<=max)))
    }
  }
  
  return(data)
}

# Check if a value (or vector of values) belongs in an interval
between <- function(evi, intervals){
  
  lower = intervals[1]
  upper = intervals[2]
  
  if(length(evi)>1){
    ans = c()
    for(i in 1:length(evi)){
      x = evi[i]
      ans[i] = between(x, intervals)
    }
    return(ans)
  }
  
  x = evi
  if((lower <= x | abs(lower-x)<10^-16) & (x <= upper | abs(upper-x)<10^-16)){
    ans = TRUE
  }else{
    ans = FALSE
  }
  return(ans)
}

# check if a specific density is valid for a particular evidence (or set of evidences)
keep_fx_vectorize <- function(bn, evi, parentIntervals){
  # bn <- bayesian network
  # evi <- evidencias
  # parentIntervals <- intervalos de evi en un data.frame
  
  # res_fin <- c()
  # Almacenamos la cantidad de intervalos
  nInt = nrow(parentIntervals)
  
  # Inicializar data.frame para cada intervalo
  res_fin <- data.frame(matrix(nrow = nrow(evi), ncol = nInt, dimnames = list(NULL, paste0('Int_',1:nInt))))
  # Almacenamos los padres que presentan evidencia
  n = colnames(evi)[which(colnames(evi)%in%colnames(parentIntervals))]
  
  # Se almacenan los niveles (o dominios) de los padres que presentan evidencias
  st = levels(bn)[n]
  
  special = c()
  single_cases = c()
  for(i in n){
    # busca los padres cuyos intervalos sean un solo punto
    equal.bounds = sapply(parentIntervals[,i], function(x){length(unique(x))==1 & is.numeric(x)})
    if(any(equal.bounds)){
      special = c(special, i)
      # Determine which intervals are left-closed
      # other parents intervals
      opi = parentIntervals[which(names(parentIntervals)!=i)]
      # opic = sapply(opi, as.character)
      opic = sapply(1:nInt, function(i){paste(unlist(opi[i,]), collapse = '')})
      opi_unique = unique(opic)
      opic_dup = match(opic, opi_unique)
      freq_opic = table(opic_dup)
      # single cases
      single_opic = which(opic %in% opi_unique[which(freq_opic==1)])
      # single_cases = unique(c(single_cases, single_opic))
      single_cases[[i]] = single_opic
    }
  }
  
  for(k in 1:nInt){
    # Se queda con el intervalo k
    obs = parentIntervals[k,, drop = FALSE]
    # Inicializa para almacenar el resultado
    res <- data.frame(matrix(nrow = nrow(evi), ncol = length(n), dimnames = list(NULL, n)))
    # Buscar el intervalo donde se situa la observacion
    for(i in 1:length(n)){
      # Almacenar los intervalos del padre n[i]
      intervals = unlist(obs[,n[i]])
      if(is.numeric(intervals)){
        # Almacenamos el estado minimo de la variable n[i]
        minx = min(st[[n[i]]])
        maxx = max(st[[n[i]]])
        # Almacenamos la evidencia de n[i]
        x = evi[,n[i]]
        # Minimo y maximo del intervalo del padre n[i]
        lower = intervals[1]
        upper = intervals[2]
        # add code for special case of zero-inflation
        if(n[i] %in% special){
          criteria =  ((lower==upper) & (lower == minx))|(k%in%single_cases[[n[i]]])
        }else{
          # ¿Es el extremo inferior de la particion el mismo valor que el minimo de la variable?
          criteria = abs(lower-minx)<10^-12
          # Modificar el extremo superior si se corresponde con el valor maximo de la variable
          if(abs(upper-maxx)<10^-12){
            upper = upper + 0.0001 
          }
        }
        
        if(criteria){
          # 'lower' is the lower bound of the interval and matches the smallest value observed in the data
          # Sometimes, R returns a difference of <10^-13, which results in a false "ans = FALSE"
          # To fix it, lower is replaced by lower-0.0001
          lower = lower-0.0001
          ans = lower<=x & x<=upper
        }else{
          ans = lower<x & x<=upper
        }
      }else{
        x = evi[,n[i]]
        ans = x == intervals
      }

      
      res[i] <- ans
    }
    res_fin[k] <- rowSums(res)==i
    # res_fin[k] <- all(res)
  }
  
  return(res_fin)
}

# Simplify parent domain of a single node.
# For instance, a node with 1 parent has 3 cases:
# 1 valid for the entire domain of the parent and the other 2 
# for 2 different intervals. The function returns which cases can be 
# summed and its domain.
simplifyDomain = function(dom){
  vars = colnames(dom)
  
  A = list()
  for(r in 1:nrow(dom)){
    df = data.frame(matrix(ncol = length(vars), nrow = 2, dimnames = list(NULL, vars)))
    for(c in 1:ncol(dom)){
      df[,c] = dom[[r,c]]
    }
    df = unique(df)
    A[[r]] = df
  }
  A
  
  if(is.numeric(dom)){
    # check which MOPs can be summed (those defined for same intervals)
    # 1. Get interval break points from the domain 
    breaks = sort(unique(unlist(dom)))
    
    # 2. Create intervals from breaks
    intervals = lapply(1:(length(breaks)-1), function(i){
      x = as.data.frame(breaks[i:(i+1)])
      colnames(x) = vars
      x
    })
    
    # 3. For each case, find the interval(s) where it overlaps
    tosum = list()
    
    for(i in 1:length(intervals)){
      tosum[[i]] = which(unlist(sapply(A, FUN = intersect2matrices, intervals[[i]])[1,]))
    }
    newDom = intervals 
  }else{
    #  check which MOPs can be summed (those defined for same intervals)
    tocheck = list()
    for(i in 1:length(A)){
      tocheck[[i]] = which(unlist(sapply(A, FUN = intersect2matrices, A[[i]])[1,]))
    }

    # remove duplicates
    tocheck = unique(tocheck)

    # for each case, check if all the domains overlap
    dom.intersec = lapply(tocheck, function(x){intersectMatrices(A[x])})

    # If TRUE, the domains overlap
    tosum = tocheck[which(sapply(dom.intersec, '[[', 'logical')==TRUE)]

    newDom = lapply(dom.intersec[which(sapply(dom.intersec, '[[', 'logical')==TRUE)], '[[', 'intersection')
  }
  
  return(list(tosum = tosum, newDom = newDom))
}



# Match several node cases according to the parent (target) domain.
# The function returns which cases can be multiplied and its domain.
matchTargetSplits <- function(splits){
  
  n.intervals = sapply(splits, nrow)
  
  # posibles combinaciones de los intervalos de cada nodo
  intervals = lapply(n.intervals, function(x){seq(1:x)})
  combination.intervals = expand.grid(intervals)
  
  # comprobar si cada combinacion de intervalos se puede multiplicar
  # y obtener el dominio sobre el cual la funcion es valida
  
  functions = splits
  
  k = 1
  selected = list()
  res = list()
  
  for(i in 1:nrow(combination.intervals)){
    comb = combination.intervals[i,, drop = FALSE]
    
    for(j in 1:length(comb)){
      selected[[j]] = functions[[j]][unlist(comb[j]),, drop = FALSE] # intervalo dentro del nodo
      selected[[j]] = cbind(selected[[j]], 'x')
      
    }
    selected
    
    # combined factor's domain
    domain = getFactorDomain(selected); domain
    
    if(!is.null(domain)){
      
      CPDs = comb
      attr(CPDs, 'DomainFactorParents') = domain
      res[[k]] = CPDs
      k = k+1
      
    }
    
  }
  return(res)
}

# Evaluate the BN on a data.frame of evidences.
# It works if all variables but one are observed.
# The function returns a data.frame containing the product of
# all the variables for each evidence
setEvidence_vectorized <- function(bn, target, evidence){
  # bn es la red bayesiana
  # target son las variables a predecir (string)
  # evidence son las observaciones con las que se predicen sin los target
  
  # browser()
  # Calcula la cantidad de valores a predecir
  N = nrow(evidence)
  # Calcula las variables observadas
  obs = names(bn)[names(bn)%in%names(evidence)]
  # target = names(bn)[!(names(bn)%in%names(evidence))]
  
  # Determina cuales son los padres de las variables Target
  logiTargetChild = sapply(lapply(bn, '[[', 'parents'), # Determina los padres
                           function(x){target%in%x} # Comprueba si es padre de las variables observadas
                           )
  # Se queda con los hijos que estan observados
  logiTargetChild = logiTargetChild[obs]
  
  # # variables involved in each MOP (node and parents). Scope of each factor
  # factorScope = lapply(bn, attr, "scope")
  # parents = lapply(factorScope, function(x){x[-length(x)]})
  # 
  # 
  # 
  
  eval_functions = list()
  
  
  allObsNode = c() # index of node with all parents observed
  targetChild = c() # index of node with target as parent
  
  # Bucle que pasa por cada variable observada
  for(i in 1:length(obs)){
    # Almacenamos la distribucion de probabilidad de la variable observada
    observedNode = bn[[obs[i]]]
    # Almacenamos los padres de la observacion 
    parents = observedNode$parents
    
    # almacenamos la funcion de densidad en fx en un data.frame
    if(observedNode$type == 'Continuous'){
      fx = observedNode$functions
    }else{
      fx = collapseDiscreteCPD(observedNode)
    }
    fx
    # encontrar densidades validas para valor observado de padres -----#
    # te quedas con los padres que son obsevarvados 
    # Pregunta: ¿Por que no se hace justo despues de la linea 274?
    parents = parents[parents%in%obs]
  # browser()
    # Pregunta: ¿if(length(parents)>0) es equivalente a que no sea NULL?
    if(length(parents)>0){ # si la cantidad de padres observados es mayor que 0
      # Se queda con las restricciones de todos los padres de la variable observada
      A = fx[observedNode$parents]
      # Almacenamos en evi la evidencia de los padres observados
      evi = evidence[parents]
      
      # (intervalos cerrados por la derecha: a<x<=b)
      # Devuelve en que dominio de los padres se situa la evidencia 
      keep = keep_fx_vectorize(bn, as.data.frame(evi), A)
      
    }else{
      keep = data.frame(matrix(data = TRUE, ncol = nrow(fx), nrow = N, dimnames = list(NULL, paste0('Int_',1:nrow(fx)))))
    }
    # ------------------------------------------------------#
    
    # evaluar mop con valor observado ----------------------#
    evi = evidence[observedNode$node]
    
    n = ncol(fx)
    
    if(observedNode$type == 'Continuous'){
      keep_functions = keep
      for(k in 1:nrow(fx)){# each row of fx corresponds to a function
        # the function is stored in the last column of object 'fx'
        evalInt = new_mop(eval.motbf(fx[k,n][[1]], evi))
        keep_functions[k] = keep[k]*evalInt
      }
    }else{# Discrete nodes
      keep_functions = keep
      
      for(k in 1:nrow(fx)){
        s = fx[[observedNode$node]][[k]];# Estados
        p = fx[[n]][[k]];# Probabilidades
        
        skeep = apply(outer(unlist(evi), s, FUN = "=="),1, which); skeep
        
        pkeep = p[skeep]
        
        keep_functions[k] = keep[k] *  p[skeep]
      }
      
    }
    # ------------------------------------------------------#
    
    # Simplificar keep_functions ---------------------------#
    # keep_functions tiene tantas columnas como filas el fx.
    # Si todos los padres están observados, cada fila solo será TRUE en una columna,
    # por lo que podemos sumar por filas y obtner una única columna: esto es
    # la densidad del nodo para cada evidencia
    
    if(logiTargetChild[i]==FALSE & all(rowSums(keep) <= 1)){# all parents are observed
      # If rowSums(keep) == 0, the combination of parents observed in the evidence did not occur in the training set
      allObsNode = c(allObsNode, i)
      keep_functions = data.frame(Int_1 = rowSums(keep_functions))
      
    }else{# one parent is the class (not observed)
      # sum columns that belong to the same class interval
      targetChild = c(targetChild, i)
      targetIntervals = fx[target]

      simpDom = simplifyDomain(targetIntervals)

      ids = simpDom$tosum
      dom = simpDom$newDom
      if(bn[[target]]$type == 'Discrete'){
        dom = lapply(dom, unique)
      }
      simpTargetIntervals = as.data.frame(cbind(dom))
      
      keep = data.frame(matrix(ncol = length(ids), nrow = N,
                               dimnames = list(NULL, paste0('Int_',1:length(ids)))))
      
      for(j in 1:length(ids)){
        id = ids[[j]]
        keep[j] = rowSums(keep_functions[id])
      }
      
      keep_functions = keep
      attr(keep_functions, 'targetIntervals') = simpTargetIntervals
      
    }
    # ------------------------------------------------------#
    
    # Guardar resultados de cada nodo
    eval_functions[[i]] = keep_functions
    
  }

  # Multiplicar densidad de nodos observados (los que no son la clase) -----#
  
  # Nodos cuyos padres estan observados (ninguno es la clase): resultado es un data.frame de una única columna
  evalAllObs = Reduce('*', eval_functions[allObsNode])
  
  if(is.null(targetChild)){# target variable is a leaf node (does not have children)
    attr(evalAllObs, 'targetDomain') = bn[[target]]$functions[target][1,,drop = FALSE]
    return(evalAllObs)
  }
  
  # find nodes whose parent is the target and multiply them -----#
  evalTargetChild = eval_functions[targetChild] # es una lista de tantos elementos como nodos cuyo (algun) padre es la clase
  targetIntervals = lapply(evalTargetChild, attr, 'targetIntervals')
  
  # if the intervals are the same, just multiply the nodes as they are
  if(all(sapply(targetIntervals, identical, targetIntervals[[1]]))){
    evalTargetChildRed = Reduce('*', evalTargetChild)
    domInter = targetIntervals[[1]]
    
  }else{# if the target intervals are different, find the domain intersection
    
    splits = matchTargetSplits(targetIntervals)
    domInter = sapply(splits, attr, 'DomainFactorParents')
    domInter = as.data.frame(cbind(domInter))
    
    evalTargetChildRed = data.frame(matrix(ncol = length(splits), nrow = N))
    
    for(k in 1:length(splits)){
      sp = splits[[k]];sp
      sel = mapply(function(x, y) x[y], evalTargetChild, sp)
      evalTargetChildRed[k] = Reduce('*', sel)
    }
  }
  
  #--------------------------------------------------------------#
  if(is.null(evalAllObs)){# if the class is parent of all nodes
    result = evalTargetChildRed
  }else{
    ints = ncol(evalTargetChildRed)
    
    temp = data.frame(rep(evalAllObs, ints))
    
    temp2 = list(evalTargetChildRed, temp)
    
    result = Reduce('*', temp2)
  }

  
  attr(result, 'targetDomain') = domInter
  return(result)
}

check.names = function(bn, evidence){
  # Check names
  if(!all(names(evidence)%in%names(bn))){
    unkownVariable = names(evidence)[which(!(names(evidence)%in%names(bn)))]
    stop("Some evidence names do not match any model names: ", paste(unkownVariable, collapse = ", "))
  }else{
    TRUE
  }
}

check.values = function(bn, evidence){
  res = c()
  
  for(i in 1:length(bn)){
    node = bn[[i]]$node
    if(node%in%names(evidence)){
      domain = levels(bn)[[node]]
      evi = evidence[[node]]
      if(bn[[i]]$type == 'Discrete'){
        if(!any(evi == domain)){
          stop("The evidence value of ", node," is not included in its set of possible values: ", paste0('{', paste(domain, collapse = ", "), '}'))
        }else{
          res[i] = TRUE
        }
      }else{
        checkEviValue = between(evi, domain)
        if(!all(checkEviValue)){
          stop("The evidence value ", paste(evi[which(checkEviValue == FALSE)], collapse = ', ')," of ", node," is outside its domain: ", paste0('[', domain[1], ', ', domain[2], ']'))
        }else{
          res[i] = TRUE
        }
      }
      
    }
  }
  return(all(res, na.rm = TRUE))
}

.pred = function(eval, fx){
  
  
  # browser()
  domainTarget = as.list(attr(eval, 'targetDomain'))[[1]]
  n = ncol(eval)
  N = nrow(eval)
  
  
  # if(is.null(bn[[target]]$children)){
  #   # If target is a leaf node, and since 'data' is filtered for a combination of parents, 
  #   # we can just compute the expected value of the target density for that combination of parents. 
  #   # evidf <- data[1, which(colnames(data)!=target), drop = FALSE]
  #   # bnq = setEvidence_VE(bn, as.list(evidf))
  #   # fx_node = bnq[[target]]$functions
  #   # fx_target = fx_node[,ncol(fx_node)][[1]]
  #   # E = expectedValue(fx_target)
  #   E = expectedValue(fx)
  #   return(rep(E, N))
  #   
  # }
  
  splitFx = splitMOP(fx)
  
  coefFx = splitFx$Coefficients
  target = splitFx$Variables
  
  b = ifelse(splitFx$Exponents == "", 0, splitFx$Exponents)# exponente para término independiente
  b = gsub(paste0('\\*',target, '|\\^'),'',b)
  b = as.numeric(ifelse(b == '', 1, b))
  
  # mutiplicar todas las observaciones por densidad target
  Q = list()
  ak = list()
 
  for(k in 1:n){
    # p = sapply(eval, '[[',k)
    p = eval[,k]
    dom = unlist(domainTarget[[k]])
    
    mcoef = outer(p, coefFx)
    Q[[k]] = mcoef
    
    # integrar (para normalizar)
    w = apply(mcoef, 1, function(i){i/(b+1)}) # coeficientes
    w = t(w)
    
    # evaluar el polinomio en el minimo y maximo 
    W1 = w%*%dom[1]^(b+1)
    W2 = w%*%dom[2]^(b+1)
    
    # sumar por filas: constantes de normalizacion
    ak[[k]] = W2-W1
  }# end loop
  
  # constante de normalizacion de cada observacion
  a = Reduce('+', ak)
  
  
  W = list()
  
  for(k in 1:n){
    # normalizar coeficientes
    P = sapply(1:N, function(i){Q[[k]][i,]/a[i]})
    P = t(P)
    
    dom = unlist(domainTarget[[k]])
    
    # Calcular esperanza: integrar (x*fx)
    w = apply(P, 1, function(i){i/(b+2)}) # coeficientes
    w = t(w)
    
    # evaluar el polinomio en el minimo y maximo 
    W[[k]] = w%*%dom[2]^(b+2) - w%*%dom[1]^(b+2)
  }
  
  # Valor esperado
  E = Reduce('+', W)
  
  
  return(E)
}



pred_vectorize = function(bn, target, data){
  # bn is bayesian network
  # target objetive variables
  # data observations
  
  # eliminamos la observaciones de las variables objetivo
  data <-data[, which(colnames(data)!=target), drop = FALSE]
  row.names(data)=NULL
  
  # Set evidence
  eval_full = setEvidence_vectorized(bn, target, data)
  
  if(bn[[target]]$type =='Discrete'){
    # Probability mass function of discrete target node
    sts =  levels(bn)[[target]]
    fx_target = collapseDiscreteCPD(bn[[target]])
    prob = data.frame(matrix(ncol = length(sts), nrow = nrow(data),
                             dimnames = list(NULL, sts)))
    
  }else{
    # Density functions of continuous target node
    fx_target = bn[[target]]$functions
  }

  ndf = ncol(fx_target)
  
  # Subset 'data' based on the target parents values/intervals
  subsetCriteria = fx_target[,-c(ndf,ndf-1), drop = FALSE]
  
  targetParents = colnames(subsetCriteria)
  
    x = c() # save results
    for(i in 1:nrow(subsetCriteria)){
      if(length(targetParents)==0){
        dataPred = data
      }else{
        sc = subsetCriteria[i,];sc
        dataPred = filterdata(data, targetParents, sc)
      }
      
      
      id = as.numeric(rownames(dataPred))
      N = nrow(dataPred)
      
      if(N>0){
        fx = fx_target[i,ndf][[1]]
        
        if(is.null(bn[[target]]$children)){# If target node is a leaf, compute the expected value directly
          
          if(bn[[target]]$type =='Discrete'){
            names(fx) = fx_target[1,ndf-1][[1]]
            E = names(which.max(fx))
            prob[id, ] = as.list(as.vector(fx))
            
          }else{
            E = expectedValue(fx)
          }
          expVal = rep(E, N)
          
        }else{ # Target node is not a leaf
          eval = eval_full[id,, drop = FALSE]
          if(bn[[target]]$type =='Discrete'){
            pd = eval*as.list(fx)
            # normalize posterior distribution
            pdN = pd/rowSums(pd)
            prob[id,] = pdN
            expVal = sts[apply(pdN,1, which.max)]
          }else{
            expVal = .pred(eval, fx)
          }
          
        }

        x[id] = expVal
      }
      
    }


  if(bn[[target]]$type =='Discrete'){
    x = factor(x, levels = sts)
    attr(x, 'prob') = prob
  }
  return(x)
}


#' Predict from an MoTBF Bayesian Network
#' @param object An object of class "motbf_fit"
#' @param target A character string indicating the variable to be predicted
#' @param data A data frame containing the predictive variables
#' @param method A character string indicating the method to perform the inference. Options are NULL, "ve" (variable elimination) and "fs" (forward sampling).
#' @param prob A logical value. If TRUE, the probabilities used for classification are returned as an attribute.
#' @param parallel A logical value. If TRUE, parallelization is carried out.
#' @param ... Additional arguments used if method = "fs".
#' @importFrom parallel detectCores mclapply makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
predict.motbf_fit<-function(object, target, data, method = NULL, prob = FALSE, parallel = FALSE, ...){

  evidf<-data[, which(colnames(data)!=target), drop = FALSE]
  
  # if all variables but the target are observed, and argument 'method' is null
  # use the vectorized function for prediction
  if(all(names(levels(object))[which(names(levels(object))!= target)] %in%colnames(evidf)) & is.null(method)){
    result = pred_vectorize(bn = object, target = target, data = evidf)
    
    if(prob == FALSE & object[[target]]$type=='Discrete'){
      attr(result, 'prob')= NULL
    }
    
    return(result)
  }
  
  N = nrow(evidf)
  
  # States/domain of target variable
  sts<-levels(object)[[target]]
    
    # make function to call variable elimination (ve) and forward sampling (fs)
    ve <- function(i){
      a = variableElimination(object, target, evidence = as.list(evidf[i, ,drop = FALSE]))
      
      if(object[[target]]$type=='Continuous'){
        b = expectedValue(a)
      }else{
        b = as.numeric(unlist(a[,length(a)]))# extract probability distribution (it is always in the last column)
        
        names(b) = a[,target]# assign class name to each value
        
        b = b[match(sts, names(b))]# make sure the class order is correct
      }

      
      return(b)
    }
    
    
    fs <- function(i,...){
      a = get_approx_posterior(bn = object, target = target, evidence = evidf[i,], parallel = parallel,...)
          
      if(object[[target]]$type=='Continuous'){
        b = expectedValue(a$fx)
      }else{
        b = a$fx$coeff
      }
      
      return(b)
    }
    
    if(method=="fs"){
      method = fs
    }else if(method=="ve" | is.null(method)){
      method = ve
    }else{
      stop("Method not found")
    }
    
    
    if(parallel == TRUE & .Platform$OS.type == "unix"){ # Parallelize code in macOS
      
      cores = detectCores()-1
      
      res <- mclapply(1:N, function(i){
        method(i)
      }, mc.cores = cores)
      

    }else if(parallel == TRUE & .Platform$OS.type != "unix"){# Parallelize code in any other system different from macOS
      
      # setup parallel backend to use many processors
      cores = detectCores()
      cl <- makeCluster(cores[1]-1) #not to overload your computer
      registerDoParallel(cl)
      
      

      #####################################-
      # res <- foreach(i = 1:N, .packages=c('bnlearn','ggm'), .export = foo) %dopar% {
      #   method(i)
      # }
      i <- NULL
      res <- foreach(i = 1:N, .packages=c('bnlearn','ggm', 'MoTBFs')) %dopar% {
        method(i)
      }
      stopCluster(cl) #stop cluster
      
    }else{ # DO NOT parallelize code
      
      
      res <- lapply(1:N, function(i){
        method(i)
      })

    }
    
  if(object[[target]]$type=='Discrete'){ ## Predict discrete variable

    dist = as.data.frame(do.call(rbind, res))
    result = names(dist)[apply(dist, 1, which.max)]
    
    result = factor(result, levels = sts)
    
    if(prob==TRUE){
      attr(result,'prob')<-dist
    }
    
  }else{ ## Predict continuous variable
    result = unlist(res)
  }
  
  return(result)
}


