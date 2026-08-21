# if(!require(stringr))
# { install.packages("stringr")
#   library(stringr)}

# if(!require(mgsub))
# { install.packages("mgsub")
#   library(mgsub)}
# 
# 
# if(!require(ggm))
# { install.packages("ggm")
#   library(ggm)}




#* The elimination order can be manually specified, using argument elimOrder, 
#* or automatically obtained, using the topological order in the DAG. 
#* If elimOrder is not specified, the topological order is computed.
#* evidence is a data.frame or list
#* target is a string containing the name of the variable in the network
#* bn is a Bayesian network, the object returned by MoTBFs_Learning().
#* 
#* 
#* Barren node:  A barren node is a node with no children, and that is not a 
#* findings node or a target node.  During belief updating and finding optimal 
#* decisions, nodes that are barren don’t influence the results, and may simply be removed.

#' Exact inference
#' 
#' Compute the posterior distribution of a variable of interest given some evidence. 
#' The variable elimination algorithm is used.
#' 
#' @param bn An object of class \code{motbf_fit}, obtained from function \link{motbf.fit}.
#' @param target A character string equal to the name of the variable of interest.
#' @param evidence A \code{data.frame} of one row containing the value of the observed variables. A list can also be provided.
#' @param elimOrder The elimination order can be manually specified as a vector containing the 
#' names of the variables, in the desired order. If elimOrder is not specified, the topological order is computed.
#' @export
#' @return The posterior probability distribution of the target variable as an object of 
#' class \code{univmotbf} or \code{piecewisemop} if \code{target} is continuous, 
#' or a matrix if \code{target} is discrete.
#' @importFrom ggm topOrder
#' 
#' @examples
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
#'   obs <- data.frame(lip = "0.48", alm1 = 0.55, stringsAsFactors=FALSE)
#'   node <- "alm2" 
#' ve = variableElimination(bn, target = node, evidence = obs)

variableElimination <- function(bn, target, evidence = NULL, elimOrder = NULL){
  
  if(is.null(bn)){
    stop('Argument "bn" is empty with no default. An object returned by MoTBFs_Learning() must be provided')
  }
  if(!is.motbf_fit(bn)){
    stop("An object of class 'motbf_fit' must be provided")
  }
  
  # Store BN object in a convenient formatted list
  # bn2 = getFormatedBN(bn)
  
  if(is.null(target)){
    stop('A target node must be specified')
  }
  if(!any(names(bn)%in%target)){
    stop('Target name given is not in the node set. Check name.')
  }
  # browser()
  if(is.data.frame(evidence)){
    if(nrow(evidence)>1){
      warning("VariableElimination() supports one piece of evidence only. Only first row of data.frame 'evidence' will be used.
              To obtain a numerical prediction of each row of 'evidence', use method 'predict()' instead.")
      evidence = evidence[1,]
    }
    evidence = as.list(evidence)
  }else if(is.list(evidence)){
    if(any(sapply(evidence, length)>1)){
      warning("VariableElimination() supports one piece of evidence only. Only first case of 'evidence' will be used.
              To obtain a numerical prediction of each case of 'evidence', use method 'predict()' instead.")
      evidence = lapply(evidence, '[[', 1)
    }
  }
  
  # save original BN
  inputBN = bn
  
  
  # Save target type (discrete/continuous)
  discreteTarget = bn[[target]]$type == 'Discrete'
  
  # Remove barren nodes. 
  bn = removeBarrenNodes(bn, target , evidence )
  
  
  # d-separation based pruning. TO DO ----
  
  
  # variables involved in each MOP (node and parents). Scope of each factor
  factorScope = lapply(bn, attr, "scope")
  
  
  if(is.null(elimOrder)){
    elimOrder = getTopoOrder(factorScope)
    elimOrder = elimOrder[length(elimOrder):1]
  }
  
  
  
  # remove target from elimination order
  targetPos = which(elimOrder%in%target)
  if(length(targetPos)>0){
    elimOrder = elimOrder[-targetPos]
  }
  
  
  #* 1. Coger todos los MOPs de la red que contengan alguna de las variables 
  #* de la evidencia, y sustituir la variable por el valor observado.
  
  if(!is.null(evidence)){
    if(target%in%names(evidence)){
      stop('The target variable cannot be observed')
    }
    check.names(bn, evidence)
    check.values(bn, evidence)
    bn = setEvidence_VE(bn = bn, evidence = evidence)
    
    # remove evidence variables from elimination order
    evidencePos = which(elimOrder%in%names(evidence))
    
    elimOrder = elimOrder[-evidencePos]
  }

  
  
  ##********* COMIENZO BUCLE ***********
  if(length(elimOrder)>0){
    
    for(i in 1:length(elimOrder)){
      
      #  2. Para cada variable X distinta de Y, coger todos los 
      #  MOPs que contengan a X, multiplicarlos y luego borrar X marginalizando.
      X = elimOrder[i]
      
      
      # cat('****** Eliminating variable', X, '******* \n')
      # Actualizar variables involucradas en cada MOP (nodo y sus padres)
      factorScope = lapply(bn, attr, "scope")
      
      # indice de las mops que hay que multiplicar (donde aparece X)
      mop_ind = which(detect_string(factorScope, X)==TRUE)
      
      subBN = subsetNetwork(bn, mop_ind)
      
      
      # cat('****** ··· Multiplying MOPs ******* \n') # multiply -----
      
      mjoint = multiplyMOPs(subBN)
    
      
      # eliminar padres observados de los atributos del factor?
      # attr(mjoint, 'parents') = attr(mjoint, 'parents')[-which(attr(mjoint, 'parents')%in%names(evidence))]
      
      # cat('****** ··· Marginalizing', X, ' ******* \n') # marginalyze ----
      # Excluir a X de la lista de variables de mjoint y simplificar despues de marginalizar X
      
      marg1 = marginalizeFactor(mjoint, X)

      
      # format the new factor 'mm' as the original nodes
      newFactor = marg1
      newFactor$node = paste0('f',i)
      
      
      
      #  3. Sustituimos todos los MOPs que hemos combinado 
      #  en el paso 2 por el que queda después de borrar X.
      bn = replaceNode(bn, mop_ind, newFactor)

    }
    #******** FIN BLUCLE *************
  }
  
  #*5. Combinar los MOPs que queden y normalizar el resultado, 
  #*que sería la distribución sobre Y

# combinar mops -----

  if(discreteTarget){
    # functions = lapply(bn, '[[', 'functions')
    
    # Unir objetos que contienen las distribuciones
    # functions1
    # functions = functions1
    # 
    # functions = lapply(functions, function(res){
    #   nc = ncol(res)
    #   numNode = which(sapply(res[,-nc, drop = FALSE], is.list))
    #   if(length(numNode)!=0){
    #     resN = res[,numNode, drop = FALSE]
    # 
    #     res[,numNode] = sapply(resN, paste)
    #     res
    #   }else{
    #     res
    #   }
    # 
    # })
    # browser()
    # Multiply discrete nodes
    # Multiply all but one to avoid problems with domain of continuous variables
    discNodes = which(sapply(bn, '[[', 'type')=='Discrete')
    if(length(discNodes)>0){
      discNodes = discNodes[-c(1)]
    }
    contNodes = which(!(1:length(bn)%in%discNodes))
    if(length(discNodes)>0){
      bnDisc = subsetNetwork(bn, discNodes)
      
      
      functions = lapply(bnDisc, '[[', 'functions')
      fooDisc = functions[[1]]
      if(length(functions)>1){
        for(i in 2:length(functions)){
          fooDisc = merge(fooDisc, functions[[i]], sort = FALSE)
        }
      }
    }
 
    
    
    # Multiply continuous nodes (multiplyMOPs() intersects the domain appropriately)
    # contNodes = which(sapply(bn, '[[', 'type')!='Discrete')
    
    if(length(contNodes)>0){
      bnCont = subsetNetwork(bn, contNodes)
      bnCont = multiplyMOPs(bnCont)
      fooCont = bnCont$functions
    }

    
    # Merge discrete and continuous
    if(exists('fooCont') & exists('fooDisc')){
      res = merge(fooDisc, fooCont, all = TRUE)
    }else if(exists('fooCont')){
      res = fooCont
    }else{
      res = fooDisc
    }
    

# Extraer columnas que contienen las distribuciones
    res2 = res[,grep('CPD|Join', colnames(res)), drop = FALSE]
    
    # Convertir objeto a data.frame para poder calcular el producto por filas
    res2b = list()
    for(i in 1:ncol(res2)){
    
      A = res2[[i]]
      if(is.list(A)){
        
        res2b[[i]] = as.numeric(sapply(A, '[[',1))
        # res2b[[i]] = as.numeric(sapply(A, '[[', 'Function'))
        
      }else{
        res2b[[i]] = A
      }
    }
    
    names(res2b) = names(res2)
    res2b = data.frame(res2b)
    
    # Producto por filas
    # res3 = apply(res2b, 1, prod, na.rm = TRUE) # na.rm = T NO esta bien hecho
    
    # If any joint == NA, that combination does not appear in the BN object (nor in the data),
    # therefore, replace it by 0.
    res3 = apply(res2b, 1, prod, na.rm = FALSE)
    res3 = ifelse(is.na(res3), 0, res3)
    
    f2 = cbind(res[target], CPD = res3/sum(res3, na.rm = TRUE))
  }else{
    # browser()
    h = multiplyMOPs(bn)
    
    if(all(sapply(lapply(h$functions[names(evidence)], unique), length)==1)){
      h$functions =  h$functions[which(!(colnames(h$functions) %in% names(evidence)))]
    }
    
    
    fx = simplifyMOPs(h)
    if(ncol(fx$functions)!=2){
      # remove columns of observed nodes
      # (all rows are compatible with the evidence)
      fx$functions = fx$functions[which(!(colnames(fx$functions) %in% names(evidence)))]
      if(ncol(fx$functions)!=2){
        stop('ncol is not equal to 2')
      }
    }
    
    
    # Before normalizing, keep the domain of the target variable only (in the mop object, "Join")
    for(i in 1:length(fx$functions$Join)){
      # fx$functions$Join[[i]]$Domain = fx$functions$Join[[i]]$Domain[,target, drop = FALSE]
      # fx$functions$Join[[i]]$Domain = fx$functions[i,target, drop = FALSE]
      fx$functions$Join[[i]]$Domain = matrix(unlist(fx$functions[i,target, drop = FALSE]), dimnames = list(NULL, target))
    }

    f2=normalizeMOP(fx)
    
    if(length(f2)>1){
      f2 = new_piecewisemop(f2)
    }
    if(length(f2)==1){
      f2 = unlist(f2, recursive = FALSE)
      f2 = new_mop(f2)
    }
  }

  return(f2)
}




# FUNCIONES PARA MULTIPLICAR MOPS -----------------------------------------
# NO EXPORTADAS
# splitMOP() moved to file mop.R
# splitMOP <- function(mop){
#   # library(stringr)
#   
#   # Extraer bases y exponentes
#   str = as.character(abs(coeffMOP(mop)))
#   
#   test = tryCatch({as.character(mop$Function)}, error=function(e){NULL})
#   if(is.null(test) & is.numeric(mop) & length(mop)==1){
#     test = mop
#   }
#   # test = as.character(mop$Function)
#   
#   # for(i in 1:length(str)){
#   #   test = str_remove(test, as.character(str)[i])
#   # }
#   
#   
#   # remove sci notation
#   test = gsub('e-','', test)
#   test = gsub('e\\+','', test)
#   
#   # Split each term by + or -
#   terms = strsplit(test, split = c("(-)|(\\+)"))[[1]]
#   # If first coef is negative, remove first empty element
#   if(coeffMOP(mop)[1]<0){
#     terms = terms[-1]
#   }
#   # split by * within each term
#   terms = strsplit(terms, split = '(?<=.)(?=\\*)', perl = TRUE);terms
#   
#   # keep exponents (second element)
#   # exponents = sapply(terms, '[',2);exponents
#   # exponents = ifelse(is.na(exponents), '', exponents);exponents
#   
#   exponents = sapply(terms, grep, pattern = "\\*", value = TRUE);exponents
#   exponents = sort(unique(unlist(exponents)));exponents
#   
#   # b = strsplit(test, split = c("(-)|(\\+)"))[[1]]
#   # # in case of sci-notation number format, or coefficients in the literal part, 
#   # # find and remove them
#   # # suppressWarnings({b[!is.na(as.numeric(b))] = ""}) 
#   # d = suppressWarnings({sapply(strsplit(b, '\\*'), function(x)ifelse(!is.na(as.numeric(x)),"", x))})
#   # b = sapply(d, paste, collapse = '*')
#   # 
#   # b = c("",b[b!=""])
#   
#   return(list(Coefficients =coeffMOP(mop), Exponents = exponents )) 
# }

reduceTerms <- function(x){
  if(x == ""){
    r = c("")
    return(r)
  }
  # separar variables
  n = strsplit(x, split = "*", fixed = TRUE)[[1]][-1]
  
  check = c()
  r = c()
  for(i in 1:length(n)){
    # comprobar si cada variable aparece varias veces
    if(i %in% check){
      next
    }
    # base
    b = strsplit(n[i], split = "^", fixed = TRUE)[[1]][1]
    v = which(grepl(b, n, fixed = TRUE))
    check = c(check,v)
    p = n[v]
    if(length(p)==1){
      r[i]=p
      next
    }
    # si una variable aparece varias veces, sumar exponentes
    q = strsplit(p, split = "^", fixed = TRUE)
    
    # sumar exponentes
    grado = as.numeric(sapply(q, "[", 2))
    
    # si es grado 1, aparece como NA. Cambiar a 1
    grado = ifelse(is.na(grado),1, grado)
    addExp = sum(grado)
    
    base = q[[1]][1]
    r[i] <- paste0(base,"^",addExp)
    
  }
  
  r = as.vector(stats::na.omit(r))
  r = paste0("*",paste0(r, collapse = "*"), collapse = "")
  return(r)
}

simplifyPolynomial <- function(cof, m){
  check <- c()
  res <- c()
  signo <- c()
  for(i in 1:length(m)){
    ind <- c(i)
    if(i %in%check){
      next
    }
    for(j in 1:length(m)){
      if(i ==j){
        next
      }
      
      p = strsplit(m[i], "*", fixed = TRUE)[[1]][-1]
      
      q = strsplit(m[j], "*", fixed = TRUE)[[1]][-1]
      
      if(all(p%in%q)& all(q%in%p)){
        check <- c(check,j)
        # terminos que se pueden sumar
        ind = c(ind,j)
      }
      
    }
    # coeficientes sumados
    res <- c(res,paste0(sum(cof[c(ind)]),m[i]))
    signo <- c(signo,sign(sum(cof[c(ind)])))
  }
  sol = paste0(ifelse(signo[-1]>=0,"+", ""), res[-1], collapse = "")
  sol = paste0(res[1],sol)
  # sol = paste0(ifelse(signo>=0,"+", ""), res, collapse = "")
  return(sol)
}

multiply2Polynomials = function(fx1, fx2){
  
  # Obtener coeficientes por un lado y variables por otro
  f1 = splitMOP(fx1)
  f2 = splitMOP(fx2)
  
  # Multiplicar terminos
  a = outer(f1$Exponents, f2$Exponents, FUN = "paste0")
  
  b = outer(f1$Coefficients, f2$Coefficients, FUN = "*")
  a;b
  
  # simplificar bases
  m <- sapply(a, reduceTerms)
  names(m) <- NULL
  m
  
  # simplificar mop
  cof = as.vector(b)
  
  poly = simplifyPolynomial(cof, m)
  
  
  # rango de las variables
  # rango = cbind(fx1$Domain,fx2$Domain)
  # rango = as.data.frame(list(fx1$Domain,fx2$Domain))

  dom_fx1 = tryCatch({fx1$Domain}, error=function(e){NULL})
  dom_fx2 = tryCatch({fx2$Domain}, error=function(e){NULL})
  if(length(dom_fx1)==0){
    dom_fx1 = NULL
  }
  if(length(dom_fx2)==0){
    dom_fx2 = NULL
  }
  dom = list(dom_fx1, dom_fx2)
  # dom = list(fx1$Domain,fx2$Domain)
  dom = dom[which(!sapply(dom, is.null))]
  rango = as.data.frame(dom, check.names = FALSE)
  
  rango = rango[,unique(colnames(rango)), drop = FALSE]
  # attr(rango, "modelVars") <- unique(colnames(rango))
  
  t5 <- list(Function=noquote(poly),Domain= rango,Subclass="mop")
  nvars = unique(c(f1$Variables, f2$Variables))
  nvars = nvars[nvars!=""]

  if(length(nvars)>1){
    output <- new_jointmotbf(t5)
  }else{
    output <- new_mop(t5)
  }
  
  
  return (output)
}

multiplyPolynomials <- function(mops){
  mm <- c()
  for(i in 2:length(mops)){ 
    if(length(mm) == 0){
      mm = multiply2Polynomials(mops[[1]], mops[[2]])
    }else{
      mm = multiply2Polynomials(mm, mops[[i]])
    }
  }
  
  attr(mm, 'DomainFactorParents') = attr(mops, 'DomainFactorParents')
  # dom_parents = attr(mm, 'DomainFactorParents')
  # # dom_parents = dom_parents[,match(colnames(mm$Domain), colnames(dom_parents))]
  # 
  # names_domain = colnames(mm$Domain)
  # names_parents = colnames(dom_parents)
  # # names_domain = colnames(mm$Domain)
  # # names_parents = colnames(mops$DomainFactorParents)
  # 
  # # id = which(names_domain %in% names_parents )
  # 
  # # for(i in 1:length(id)){
  # #   id_p = which(names_parents == names_domain[id[i]])
  # #
  # #   mm$Domain[,id[i]] = mops$DomainFactorParents[[id_p]]
  # # }
  # id = which(!(names_domain %in% names_parents ))
  # dom = cbind(mm$Domain[,id, drop = FALSE], dom_parents)
  dom = attr(mm, 'DomainFactorParents')
  dom = unique(dom)
  dom = dom[,sapply(dom, is.numeric), drop = FALSE]
  if(length(mm$Domain) !=0){
    mm$Domain = as.matrix(dom) # por qué tiene que ser una matriz?
  }else{
    mm$Domain = NULL
  }
  
  # attr(mm$Domain, "modelVars") <- colnames(dom)
  return(mm)
}


multiplyMOPs <- function(cond){
  
  if(length(cond) == 1){
    return(cond[[1]])
    # res = list()
    # functions = cond[[1]]$functions
    # for(i in 1:length(functions)){
    #   
    #   p_domain = functions[[i]]$parentInterval
    #   ch_domain = functions[[i]]$Fx$Domain
    #   id = which(!(colnames(ch_domain) %in% colnames(p_domain )))
    #   dom = cbind(ch_domain[,id, drop = FALSE], p_domain)
    #   res[[i]] = functions[[i]]$Fx 
    #   res[[i]]$Domain = dom
    # }
  }else{
    # find mops that can be multiplied together
    splits = matchSplits(cond)
    res = list()
    
    for(i in 1:length(splits)){
      mops = splits[[i]]
      mmops =   multiplyPolynomials(mops)
      attributes(mmops)
      
      res[[i]] = mmops 
    }
  }
  
  scope = colnames(attr(res[[1]], 'DomainFactorParents'))
  #*************************************************
  # Add factor parents and children as attributes
  p = unique(unlist(sapply(cond, '[', 'parents')));p
  ch = unique(unlist(sapply(cond, '[', 'children')));ch
  # n = names(cond);n
  # browser()
  n = getVariablesMOP(cond); n
  # n = unlist(unique(lapply(res, getMotbfVar)))
  factor_parents = p[!(p%in%n)]
  factor_children = ch[!(ch%in%n)]
  # attr(res, 'parents') = factor_parents
  # attr(res, 'children') = factor_children
  #*************************************************
  #*
  
  
  # Store conditional probability distributions in a data.frame, including parent's domain
  CPDs = data.frame(matrix(ncol = length(scope)+1, 
                           dimnames = list(NULL, c(scope, 'Join'))))
  
  for(i in 1:length(res)){
    domain = attr(res[[i]], 'DomainFactorParents')
    
    for(k in 1:ncol(domain)){
      if(is.numeric(domain[,k])){
        A = domain[,k, drop = TRUE]
        CPDs[[i,k]] = list(A)
        CPDs[[i,k]]  = A
      }else{
        CPDs[i, k] = unique(domain[,k])
      }
    }
    CPDs[[i, 'Join']] = res[i]
    CPDs[[i, 'Join']] = res[[i]]
  }
  
  
  newFactor <- list()
  newFactor$node = 'factor'
  newFactor['parents'] = list(factor_parents)
  newFactor['children'] = list(factor_children)
  newFactor$type = ifelse(all(sapply(cond, '[[', 'type') == 'Discrete'), 'Discrete', 
                          ifelse(all(sapply(cond, '[[', 'type') == 'Continuous'), 'Continuous', 'Hybrid'))
  newFactor$subclass = 'mop'
  
  newFactor$functions  = CPDs
  
  attr(newFactor, 'scope') = scope
  # newFactor = new_motbf_fit_node(newFactor)
  return(newFactor)
}
# multiplyMOPs <- function(cond){
#   
#   if(length(cond) == 1){
#     # res = sapply(sapply(cond, '[[','functions'),'[','Fx')
#     # attr(res, 'scope') = attr(cond[[1]], 'scope')
#     # sapply(sapply(cond, '[[','functions'),'[','parentInterval')
#     # sapply(res, '[[', 'Domain')
#     res = list()
#     functions = cond[[1]]$functions
#     for(i in 1:length(functions)){
#       
#       p_domain = functions[[i]]$parentInterval
#       ch_domain = functions[[i]]$Fx$Domain
#       id = which(!(colnames(ch_domain) %in% colnames(p_domain )))
#       dom = cbind(ch_domain[,id, drop = FALSE], p_domain)
#       res[[i]] = functions[[i]]$Fx 
#       res[[i]]$Domain = dom
#     }
#   }else{
#     # find mops that can be multiplied together
#     splits = matchSplits(cond)
#     res = list()
#     for(i in 1:length(splits)){
#       mops = splits[[i]]
#       mmops =   multiplyPolynomials(mops)
#       
#       # # update domain of parent variables
#       # domain = lapply(mops, '[[',1) # parent domain in each interval of children
#       # parents = unique(unlist(sapply(domain, names)))
#       # 
# 
#       res[[i]] = mmops 
#     }
#   }
#   
#   attr(res, 'scope') = colnames(res[[1]]$Domain)
#   #*************************************************
#   # Add factor parents and children as attributes
#   p = unique(unlist(sapply(cond, '[', 'parents')));p
#   ch = unique(unlist(sapply(cond, '[', 'children')));ch
#   # n = names(cond);n
#   n = getVariablesMOP(cond); n
#   n = unique(unlist(n))
#   factor_parents = p[!(p%in%n)]
#   factor_children = ch[!(ch%in%n)]
#   attr(res, 'parents') = factor_parents
#   attr(res, 'children') = factor_children
#   #*************************************************
#   #*
#   
#   return(res)
# }

getVariablesMOP <- function(bn){
  # fx = unlist(lapply(bn, '[[', 'functions'), recursive = FALSE)
  # varsMOP = sapply(sapply(fx, '[', 'Fx'), getMotbfVar)
  # browser()
  type = sapply(bn, '[[', 'type')
  vars = names(which(type == 'Discrete'))
  
  bnC = subsetNetwork(bn, which(type != 'Discrete'))
  if(length(bnC)>0){
    # browser()
    fx = sapply(sapply(bnC, '[', 'functions'), function(x){x[,c(length(x)), drop = FALSE]})
    fx = unlist(fx, recursive = FALSE)
    varsMOP = sapply(fx, getMotbfVar)
    
    vars = c(vars, unique(varsMOP))
  }

  return(vars)
}

# marginalizeFactor <- function(jointFactor, var){
#   marg <- list()
#   for(k in 1:length(jointFactor)){
#     
#     marg[[k]] = marginal.jointmotbf(jointFactor[[k]], var)
#   }
#   
#   attr(marg, "scope") = colnames(marg[[1]]$Domain)
#   
#   #*************************************************
#   # Add factor parents and children as attributes
#   attr(marg, 'parents') = attr(jointFactor, 'parents')
#   attr(marg, 'children') = attr(jointFactor, 'children')
#   #*************************************************
#   
#   return(marg)
# }


# marginalizeFactor <- function(joinFactor, var){
#   
#   functions = joinFactor$functions
#   marg <- functions
#   
#   pos = ncol(functions)
#   for(k in 1:nrow(functions)){
#     
#     fx = functions[[k, pos]]
#     marg[[k, pos]] =  marginal.jointmotbf(fx, var)
#   }
#   
#   marg = marg[,c(which(colnames(marg)%in%var), pos)] # remove column with domain of marginalized variable
#   
#   joinFactor$functions = marg
#   attributes(joinFactor)
#   
#   attr(joinFactor, "scope") = var
#   
#   return(joinFactor)
# }


marginalizeFactor <- function(joinFactor, X){
  
  var = attr(joinFactor, "scope")[which(attr(joinFactor, "scope") != X)]
  
  functions = joinFactor$functions
  if(is.numeric(unlist(functions[,X]))){ #variable to be marginalized is numeric
    marg <- functions
    
    pos = ncol(functions)
    for(k in 1:nrow(functions)){
      
      fx = functions[[k, pos]]
      marg[[k, pos]] =  marginal.jointmotbf(fx, var)
    }
    
    marg = marg[,c(which(colnames(marg)%in%var), pos)] # remove column with domain of marginalized variable
    
    joinFactor$functions = marg
    attributes(joinFactor)
    
    attr(joinFactor, "scope") = var
    
    mm = simplifyMOPs(joinFactor)
    
  }else{# variable to be marginalized is discrete
    f_removed = functions[,-which(colnames(functions)==X)]
    joinFactor$functions = f_removed
    mm = simplifyMOPs(joinFactor)
    
    attr(mm, "scope") = var
    
    mm$parents = setdiff(mm$parents, X)
    mm$children = setdiff(mm$children, X)
    
    
  }
  
  joinFactor = mm
  
  return(joinFactor)
}

simplify2MOPs <- function(fx1, fx2){
  #sumar mops despues de marginalizar
  mop1 = splitMOP(fx1)
  mop2 = splitMOP(fx2)
  
  # obtener terminos que se pueden sumar 
  ind1 = which(mop1$Exponents%in%mop2$Exponents)
  ind2 = which(mop2$Exponents%in%mop1$Exponents)
  if(length(ind1)!=length(ind2)){
    stop("Different length")
  }
  # ordenar 1 de los vectores
  # ind1 = match(mop1$Exponents,mop2$Exponents)
  # mop1$Exponents[ind1]
  # mop2$Exponents[ind2]
  # match(mop1$Exponents[ind1],mop2$Exponents[ind2])
  
  # ind1 = match(mop2$Exponents[ind2],mop1$Exponents[ind1])
  
  mop1exp = unlist(mop1$Exponents[ind1])
  mop1coef = mop1$Coefficients[ind1]
  
  newOrder = match(mop2$Exponents[ind2],mop1$Exponents[ind1])
  mop1exp <- mop1exp[newOrder]
  mop1coef <- mop1coef[newOrder]
  mop2exp = unlist(mop2$Exponents[ind2])
  if(any(mop1exp!=mop2exp)){
    stop('MOP exponents do not match')
  }
  
  # los terminos en la posicion ind1 e ind2 estan definidos sobre el mismo producto de variables
  # sumCoef = mop1$Coefficients[ind1]+mop2$Coefficients[ind2]
  # sumTerms = paste0(ifelse(sumCoef>=0,"+", ""),sumCoef,mop1$Exponents[ind1], collapse = "")
  
  sumCoef = mop1coef+mop2$Coefficients[ind2]
  sumTerms = paste0(ifelse(sumCoef>=0,"+", ""),sumCoef,mop1exp, collapse = "")
  
  sumTerms = paste0(ifelse(sumCoef[-1]>=0,"+", ""),sumCoef[-1],mop1exp[-1], collapse = "")
  sumTerms = paste0(sumCoef[1], mop1exp[1], sumTerms)
  
  # resto de terminos (los que no se suman)
  no_ind1 = which(!mop1$Exponents%in%mop2$Exponents)
  no_ind2 = which(!mop2$Exponents%in%mop1$Exponents)
  
  terms1 = paste0(ifelse(mop1$Coefficients[no_ind1]>=0,"+", ""), mop1$Coefficients[no_ind1], mop1$Exponents[no_ind1], collapse = "")
  terms2 = paste0(ifelse(mop2$Coefficients[no_ind2]>=0,"+", ""),mop2$Coefficients[no_ind2], mop2$Exponents[no_ind2], collapse = "")
  
  res = paste0(sumTerms, terms1, terms2, collapse = "")
  
  # # rango de las variables
  # rango = cbind(fx1$Domain,fx2$Domain)
  # rango= rango[,unique(colnames(rango)), drop= FALSE]
  # attr(rango, "scope") <- unique(colnames(rango))
  
  # dominio variables
  dom_fx1 = fx1$Domain
  dom_fx2 = fx2$Domain
  if(all(sapply(list(dom_fx1, dom_fx2), is.null)) | all(sapply(list(dom_fx1, dom_fx2), length)==0)){
    rango = NULL
  }else{
    dom = intersect2matrices(fx1$Domain, fx2$Domain)
    rango = dom$intersection
  }

  # browser()
  t5 <- list(Function=noquote(res),Domain = rango, Subclass="mop")
  t5 <- new_motbf(t5)
  nVar = getMotbfVar(t5)
  # output <- jointmotbf(t5)
  if(length(nVar)>1){
    output <- new_jointmotbf(t5)
  }else{
    output <- new_mop(t5)
  }
  return(output)
}

# simplifyMOPs <- function(allmops){
#   # get domain of each mop piece
#   dom = sapply(allmops, '[[', 'Domain', simplify = FALSE)
#   
#   
#   #  check which MOPs can be summed (those defined for same intervals)
#   intersect2matrices(dom[[1]], dom[[1]])
#   tocheck = list()
#   for(i in 1:length(dom)){
#     tocheck[[i]] = which(unlist(sapply(dom, FUN = intersect2matrices, dom[[i]])[1,]))
#   }
#   # remove duplicates
#   tocheck = unique(tocheck)
#   
#   # for each case, check if all the domains overlap
#   dom.intersec = lapply(tocheck, function(x){intersectMatrices(dom[x])})
#   
#   # If TRUE, the domains overlap
#   tosum = tocheck[which(sapply(dom.intersec, '[[', 'logical')==TRUE)]
#   
#   
#   # loop over 'tosum'. Sum the mops matching the index at each element of the list.
#   # If the element contains only one index, then return the mop piece as it is.
#   k = 1
#   res = list()
#   for(k in 1:length(tosum)){
#     mops <- allmops[tosum[[k]]]
#     if(length(mops)>1){
#       mm <- c()
#       for(i in 2:length(mops)){
#         if(length(mm) == 0){
#           mm = simplify2MOPs(mops[[1]], mops[[2]])
#         }else{
#           mm = simplify2MOPs(mm, mops[[i]])
#         }
#       }
#     }else{
#       mm <- mops[[1]]
#       n = tryCatch({getMotbfVar(mm)}, error=function(e){NULL})
#       if(!is.null(n)){
#         class(mm) = ifelse(length(n)== 1, 'motbf', 'jointmotbf')
#       }
#     }
#     res[[k]] <- mm
#   }
#   attr(res, "scope") = colnames(res[[1]]$Domain)
#   
#   #*************************************************
#   # Add factor parents and children as attributes
#   attr(res, 'parents') = attr(allmops, 'parents')
#   attr(res, 'children') = attr(allmops, 'children')
#   #*************************************************
#   
#   return(res)
# }


simplifyMOPs <- function(allmops){
  # get domain of each mop piece

  dom = allmops$functions[,c(1: ncol(allmops$functions)-1), drop = FALSE]

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
  
  # Store conditional probability distributions in a data.frame, including parent's domain
  CPDs = data.frame(matrix(ncol = length(vars)+1, 
                           dimnames = list(NULL, c(vars, 'Join'))))
  
  
  # loop over 'tosum'. Sum the mops matching the index at each element of the list.
  # If the element contains only one index, then return the mop piece as it is.


  functions = allmops$functions
  
  pos = ncol(functions)
  res = functions
  fx = functions[,pos, drop = FALSE]
  if(length(tosum)>0){
    for(k in 1:length(tosum)){
      mops <- fx[tosum[[k]],]
      if(length(mops)>1){
        mm <- c()
        for(i in 2:length(mops)){
          if(length(mm) == 0){
            mm = simplify2MOPs(mops[[1]], mops[[2]])
          }else{
            mm = simplify2MOPs(mm, mops[[i]])
          }
        }
      }else{
        mm <- mops[[1]]
        n = tryCatch({getMotbfVar(mm)}, error=function(e){NULL})
        if(!is.null(n)){
          # class(mm) = ifelse(length(n)== 1, 'motbf', 'jointmotbf')
          if(length(n)==1){
            mm <- new_mop(mm)
          }else{
            mm <- new_jointmotbf(mm)
          }
        }
      }
      
      
      # replace functions by simplified one
      # for(j in 1:length(tosum[[k]])){
      #   res[tosum[[k]][[j]],pos][[1]] = list(mm)
      # }
      
      domain = newDom[[k]]
      
      for(v in 1:ncol(domain)){
        if(is.numeric(domain[,v])){
          A = domain[,v, drop = TRUE]
          CPDs[[k,v]] = list(A)
          CPDs[[k,v]]  = A
        }else{
          CPDs[k, v] = unique(domain[,v])
        }
      }
      CPDs[[k, 'Join']] = list(mm)
      CPDs[[k, 'Join']] = mm
      
    }
    # res = unique(res)
    # 
    # allmops$functions = res
    
    allmops$functions = CPDs
  }

  return(allmops)
}




intersectMatrices <- function(set){
  if(length(set)<2){
    mm = list(logical = TRUE, intersection = set[[1]])
  }else{
    mm <- c()
    for(i in 2:length(set)){
      if(length(mm) == 0){
        mm = intersect2matrices(set[[1]], set[[2]])
      }else{
        if(is.null(mm$intersection)){
          return(list(logical = FALSE, intersection = NULL))
        }
        mm = intersect2matrices(mm$intersection, set[[i]])
      }
    }
  }
  
  return(mm)
}

intersect2matrices <- function(A, B){
  # A = as.matrix(A)
  # B = as.matrix(B)
  r = c()
  
  for(i in 1:ncol(A)){
    if(is.numeric(A[,i])){
      r[i] = intersection(A[,i], B[,i])$logical
    }else{
      r[i] = all(A[,i] == B[,i])
    }
    
  }
  
  if(all(r)){
    res = TRUE
    r.m = data.frame(matrix(nrow = 2, ncol = ncol(A), 
                            dimnames = list(NULL, colnames(A))))
    for(i in 1:ncol(A)){
      if(is.numeric(A[,i])){
        r.m[,i] = intersection(A[,i], B[,i])$intersection
        # attributes(r.m) = attributes(A)
      }else{
        r.m[,i] = A[,i]
      }
      
    }
  }else{
    res = FALSE
    r.m = NULL
    
  }
  return(list(logical = res, intersection = r.m))
}

intersection <- function(x, y, value = NULL){
  if(!is.null(value)){
    if(!(value%in%c('int', 'logi'))){
      stop('"value" can only take "int" or "logi')
    }
  }else{
    value = 'both'
  }
  a = max(c(min(x),min(y)))
  b = min(c(max(x),max(y)))
  
  if(a<b){
    res = c(a,b)
    names(res) = c('min', 'max')
    res.logi = TRUE
  }else{
    res = NULL
    res.logi = FALSE
  }
  if(value == 'int'){
    return(res)
  }else if(value == 'logi'){
    return(res.logi)
  }else{
    return(list(logical = res.logi, intersection = res))
  }
}


##*******************************************
##*# format factor to get the same structure as in the original nodes
##* mm is the output of simplifyMOPs


formatFactor <- function(mm){
  
  newFactor <- list()
  newFactor$node = 'factor'
  newFactor['parents'] = list(attr(mm, 'parents'))
  newFactor['children'] = list(attr(mm, 'children'))
  # newFactor$type = node$type # discrete, continuous, hybrid
  newFactor$subclass = mm[[1]]$Subclass
  
  
  functions <- list()
  for(i in 1:length(mm)){
    piece = mm[[i]]
    
    
    # vars.parent <- attr(mm, 'scope')[-which(attr(mm, 'scope')%in% vars.f)]
    vars.parent <-attr(mm, 'parents')
    if(length(vars.parent)==0){
      parentInterval <- NULL
    }else{
      parentInterval <- data.frame(piece$Domain[,vars.parent, drop = FALSE])
      
    }
    
    if(is.numeric(piece$Function)){
      piece$Domain <- NULL
    }else{
      vars.f <- getMotbfVar(piece)
      piece$Domain <- piece$Domain[,vars.f, drop = FALSE]
    }
    
    
    
    
    interval =  list(parentInterval = parentInterval, Fx = piece)
    functions[[i]] <- interval
  }
  
  names(functions) <- paste0('Interval_',1:length(functions))
  # functions <- list(functions = functions)
  newFactor$functions  = functions
  
  attr(newFactor, 'scope') = attr(mm, 'scope')
  
  return(newFactor)
}






#*****************************
#*****************************
# get list of splits in a sub-network
# matchSplits <- function(bn){
#   n.nodes = length(bn)
#   n.intervals = sapply(lapply(bn, '[[','functions'), length)
#   
#   # posibles combinaciones de los intervalos de cada nodo
#   intervals = lapply(n.intervals, function(x){seq(1:x)})
#   combination.intervals = expand.grid(intervals)
#   
#   # comprobar si cada combinacion de intervalos se puede multiplicar
#   # y obtener el dominio sobre el cual la funcion es valida
#   
#   functions = lapply(bn, '[[', 'functions')
#   
#   k = 1
#   selected = list()
#   res = list()
#   
#   for(i in 1:nrow(combination.intervals)){
#     comb = combination.intervals[i,, drop = FALSE]
#     for(j in 1:length(comb)){
#       node = which(names(functions) == names(comb[j])) # nodo
#       selected[[j]] = functions[[node]][unlist(comb[j])] # intervalo dentro del nodo
#       
#     }
#     sel2 = unlist(selected, recursive = FALSE)
#     domain = getDomainFactorParents(sel2)
# 
#     # Si es nulo, los intervalos de un padre continuo no intersectan o 
#     # los estados de un padre discreto son distintos. En este caso, la combinación se descarta.
#     if(!is.null(domain)){
#       sel2$DomainFactorParents = domain
#       res[[k]] = sel2
#       k = k+1
#     }
#     
#   }
#   return(res)
# }

#* return CPDs to multiply together.
#* this function checks which CPDs can be multiplied, according to their parent values
#* Argument subBN is of class "motbf_fit" "motbf.ve"
matchSplits <- function(subBN){
  
  n.nodes = length(subBN)
  n.intervals = sapply(lapply(subBN, '[[','functions'), nrow)
  
  # posibles combinaciones de los intervalos de cada nodo
  intervals = lapply(n.intervals, function(x){seq(1:x)})
  combination.intervals = expand.grid(intervals)
  
  # comprobar si cada combinacion de intervalos se puede multiplicar
  # y obtener el dominio sobre el cual la funcion es valida
  
  functions = lapply(subBN, '[[', 'functions')
  
  k = 1
  selected = list()
  res = list()
  i = 1
  j = 1
  for(i in 1:nrow(combination.intervals)){
    comb = combination.intervals[i,, drop = FALSE]
    
    
    for(j in 1:length(comb)){
      node = which(names(functions) == names(comb[j])) # nodo
      selected[[j]] = functions[[node]][unlist(comb[j]),] # intervalo dentro del nodo
      
    }
    # sel2 = unlist(selected, recursive = FALSE)
    
    # combined factor's domain
    domain = getFactorDomain(selected)
    
    if(!is.null(domain)){
      # sel2$DomainFactorParents = domain
      # res[[k]] = sel2
      
      CPDs = sapply(selected, function(x){x[[length(x)]]})
      attr(CPDs, 'DomainFactorParents') = domain
      res[[k]] = CPDs
      k = k+1
      # multiply2Polynomials(CPDs[[1]], CPDs[[2]]) 
      
      # CPDs = lapply(selected, function(x){x[[length(x)]]})
      # 
      # names(CPDs) = c('Fx', 'Fx')
      # 
      # 
      # CPDs1 = do.call(cbind, CPDs)
      # domain
      # 
      # 
      # test = lapply(domain, list)
      # test2 = do.call(cbind, test)
      # test3 = as.data.frame(cbind(test2, CPDs1))
      
    }
    
  }
  return(res)
}

getFactorDomain = function(pieces){
  # Extract the parent domain in the given splits of the child nodes
  A = lapply(pieces, function(x){x[,-c(length(x)), drop = FALSE]})
  A = A[lapply(A, is.null) == FALSE]
  p = unique(unlist(lapply(A, colnames)))

  
  ls <- list()

  # For each parent, get its domain (intersection) for the given combination of splits
  for(i in 1:length(p)){
    B = lapply(A, function(x) subset(x, select = intersect(p[i], colnames(x))) )
    
    B = B[lapply(B, length)!=0];B
    
    if(isDiscreteParent(B)){
      # Check if the splits are defined for the same state of the discrete parent
      a = unlist(sapply(B, unique))
      if(length(unique(a))==1){
        ls[[i]] = c(unique(a), unique(a)) # save it twice (so that conversion into data.frame is possible when continuous parents are present)
      }else{# return NULL if the splits are not defined for the same state of the discrete parent
        return(NULL)
      }
    }else{# continuous parent
      B2 = sapply(B, as.data.frame)
      B3 = lapply(B2, unlist)
      a = max(unlist(lapply(B3,min)))
      b = min(unlist(lapply(B3,max)))
      # a = max(unlist(lapply(B,min)))
      # b = min(unlist(lapply(B,max)))
      # check if a < b. If not, the domain of the continuous parent does not intersect
      if(a<b){
        ls[[i]] = c(a, b)
      }else{
        return(NULL)
      }
    }
  }
  names(ls) = p
  domainFactorParents = as.data.frame(ls)
  return(domainFactorParents)
}
isDiscreteParent = function(B){
  res = c()
  for(i in 1:length(B)){
    x = c()
    for(k in 1:ncol(B[[i]])){
      # Check if both values are the same. If not, the parent is not discrete
      if(length(unique(B[[i]][,k]))==1){
        # Check if the values are character. If not, the parent is not discrete
        x[k] = is.character(B[[i]][,k])
      }else{
        return(FALSE)
      }
    }
    res[i] = all(x)
  }
  return(all(res))
}
getDomainFactorParents = function(pieces){
  # Extract the parent domain in the given splits of the child nodes
  A = lapply(pieces, '[[', 'parentInterval')
  A = A[lapply(A, is.null) == FALSE]
  p = unique(unlist(lapply(A, colnames)))
  
  
  ls <- list()
  # For each parent, get its domain (intersection) for the given combination of splits
  for(i in 1:length(p)){
    B = lapply(A, function(x) subset(x, select = intersect(p[i], colnames(x))) )
    
    B = B[lapply(B, length)!=0];B
    
    if(isDiscreteParent(B)){
      # Check if the splits are defined for the same state of the discrete parent
      a = unlist(sapply(B, unique))
      if(length(unique(a))==1){
        ls[[i]] = c(unique(a), unique(a)) # save it twice (so that conversion into data.frame is possible when continuous parents are present)
      }else{# return NULL if the splits are not defined for the same state of the discrete parent
        return(NULL)
      }
    }else{# continuous parent
      a = max(unlist(lapply(B,min)))
      b = min(unlist(lapply(B,max)))
      # check if a < b. If not, the domain of the continuous parent does not intersect
      if(a<b){
        ls[[i]] = c(a, b)
      }else{
        return(NULL)
      }
    }
  }
  names(ls) = p
  domainFactorParents = as.data.frame(ls)
  return(domainFactorParents)
}



#**************************************
# FUNCIONES PARA DISCRETAS  -------------------------------------
# .-ESTA FUNCION NO SE USA ----
factorProductDiscrete2 <- function(subBN){
  # All nodes are discrete
  all(sapply(subBN, '[[','type') == 'Discrete')
  # All parents are discrete
  all(lapply(levels(subBN), 'is.numeric')==FALSE)
  
  subFactorLevels = levels(subBN)
  # Combinations of nodes' states
  subFactorCombinations = expand.grid(subFactorLevels, stringsAsFactors = FALSE)
  
  
  # parents
  subFactorParents = lapply(lapply(subBN, attr, "scope"), function(x){x[-length(x)]})
  
  vars = colnames(subFactorCombinations)
  
  # Distribution of factors in "subBN"
  P = lapply(subBN, '[[', 'functions')
  i = 1
  subFactorProb = list()
  for(i in 1:nrow(subFactorCombinations)){ #loop over each combination of nodes' states
    comb = subFactorCombinations[i,]
    
    k = 1
    prob <- c()
    for(k in 1:length(P)){# Loop over factors in "subBN"
      case = P[[k]][-ncol(P[[k]])] # Last column is the probability distribution
      v = which(colnames(comb)%in%colnames(case)) # look at columns that match the scope of given factor
      combF = comb[v][match(colnames(case),colnames(comb[v]))] #make sure the column order is the same
      
      idP = which(apply( case, 1, function(x){all(x == combF)}))
      
      prob[k] = P[[k]][[ncol(P[[k]])]][idP]
      
    }
    prob
    subFactorProb[[i]] = prob
  }
  probs = do.call(rbind, subFactorProb)
  colnames(probs) = names(subFactorParents)
  
  subFactorCombinations$FactorProduct = apply(probs, 1, prod)
  
  attr(subFactorCombinations, 'class') =  c('data.frame', 'factorProduct.dicrete')
  
  
  # attr(subFactorCombinations, 'scope') = vars
  #*************************************************
  # Add factor parents and children as attributes
  p = unique(unlist(sapply(subBN, '[', 'parents')));p
  ch = unique(unlist(sapply(subBN, '[', 'children')));ch
  # n = names(cond);n
  
  n = names(subFactorParents);n
  
  factor_parents = p[!(p%in%n)]
  factor_children = ch[!(ch%in%n)]
  # attr(subFactorCombinations, 'parents') = factor_parents
  # attr(subFactorCombinations, 'children') = factor_children
  #*************************************************
  
  fnode <- list()
  fnode$node = 'newFactor'
  fnode['parents'] = list(factor_parents)
  fnode['children'] = list(factor_children)
  fnode$type = 'Discrete'
  fnode$subclass = 'Multinomial'
  fnode$functions  = subFactorCombinations
  attr(fnode, which = "scope") = vars
  
  class(fnode) = c(class(subBN),'factorProduct.dicrete')
  return(fnode)
}

#* fpd is a dataframe containing the factor product of a discrete sub network. 
#* The object is of class 'factorProduct.dicrete'
#* x is the variable to marginalize (erase) from the factor product
# .-ESTA FUNCION NO SE USA ----
# marginalizeDiscrete2 <- function(factor.d, x){
#   # if(!('factorProduct.dicrete'%in% class(factor.d))){
#   #   stop('Object "factor.d" must be of class "factorProduct.dicrete"')
#   # }
#   p = factor.d$parents
#   ch = factor.d$children
#   vars = attr(factor.d, 'scope')
#   # p = attr(factor.d, 'parents')
#   # ch = attr(factor.d, 'children')
#   # vars = colnames(factor.d)[-ncol(factor.d)]
#   fpd = factor.d$functions
#   fpd = fpd[,-c(which(vars == x)), drop = FALSE]
#   
#   # Sum out X
#   l = ncol(fpd)
#   uniqueFinalCombs = unique(fpd[,-l, drop = FALSE]) # last column is the probability distribution
#   rownames(uniqueFinalCombs)= NULL
#   
#   
#   newFactor = uniqueFinalCombs
#   # i = 1
#   for(i in 1:nrow(uniqueFinalCombs)){
#     
#     uComb = uniqueFinalCombs[i,]
#     uComb
#     check = c()
#     for(k in 1:nrow(fpd)){
#       
#       check[k] = all(fpd[k,-l] == uComb)
#       
#     }
#     fpd[which(check),]
#     sum(fpd[which(check),l])
#     newFactor[i,'SumProduct'] = sum(fpd[which(check),l])
#   }
#   
#   factor.d$functions = newFactor
#   attributes(factor.d)
#   
#   attr(factor.d, "scope") = vars[-which(vars == x)]
#   attr(factor.d, 'class') =  c("motbf.fit", "motbf.ve" , 'sumProduct.dicrete')
#   
#   
#   
#   #*************************************************
#   # Add factor parents and children as attributes
#   # attr(newFactor, 'parents') = p
#   # attr(newFactor, 'children') = ch
#   #*************************************************
#   return(factor.d)
# }







#-------------------------------------------#
# FUNCIONES PARA FORMATO DE LA RED BAYESIANA -------------------------------------



# Subset the BN keeping the attributes
subsetNetwork <- function(bn, nodes){
  if(!is.motbf_fit(bn)){
    stop('Object "bn" must be of class "motbf_fit')
  }
  factorScope = lapply(bn, attr, "scope")
  subFactorScope = unique(unlist(factorScope[nodes]))
  subFactorLevels = levels(bn)[subFactorScope]
  
  subBN = bn[nodes]
  attr(subBN, 'levels') = subFactorLevels
  attr(subBN, 'class') = attr(bn, 'class')
  
  return(subBN)
}

addNode <- function(bn, node){
  if(!is.motbf_fit(bn)){
    stop('Object "bn" must be of class "motbf_fit')
  }

  bn[[length(bn)+1]] <- node
  names(bn)[length(bn)]<- node$node

  attributes(bn)
  return(bn)
}


replaceNode <- function(bn, remove, add){
  if(!is.motbf_fit(bn)){
    stop('Object "bn" must be of class "motbf_fit')
  }
  
  bn2 = subsetNetwork(bn, -remove)
  bn3 = addNode(bn2, add)
  
  vars = unique(unlist(sapply(bn3, attr, "scope")))
  
  attr(bn3, 'levels') = levels(bn)[vars]
  
  return(bn3)
}
removeBarrenNodes <- function(bn, target, evidence){
  bnNodes = names(bn)
  possible_barren = bnNodes[!(bnNodes%in%c(target, names(evidence)))]
  if(length(possible_barren)==0){
    return(bn)
  }
  barrenNodes = names(which(sapply(sapply(bn[possible_barren], '[[', 'children'), is.null)))
  
  if(length(barrenNodes)==0){
    return(bn)
  }
  
  message(paste0('Barren nodes (', paste(barrenNodes, collapse = ", "), ") have been removed from network"))
  
  lev = levels(bn)
  
  bn2 = bn[!(names(bn)%in%barrenNodes)]
  sapply(bn2, '[[', 'children')
  for(i in 1:length(bn2)){
    ch = which(bn2[[i]]$children%in%barrenNodes)
    if(length(ch)>0){
      children = bn2[[i]]$children[-ch]
      if(length(children)==0){
        children = NULL
      }
      bn2[[i]]['children'] = list(children)
    }
  }
  sapply(bn2, '[[', 'children')
  attr(bn2, "levels") = lev[!(names(lev)%in%barrenNodes)]
  
  attr(bn2, "class") = class(bn)
  
  # recursively, remove barren nodes
  if(length(barrenNodes)>0){
    bn2 = removeBarrenNodes(bn2, target, evidence)
  }
  
  
  return(bn2)
}


#**************************************
# FUNCIONES PARA VARIABLE ELIMINATION -------------------------------------
#* Return the BN once the variables have been substituted by the observed values 
#* Works for continuous BNs only
# .-ESTA FUNCION NO SE USA ----
setEvidence <- function(bn, evidence){
  
  #*** COMPROBAR EVIDENCIA ******
  # Check names
  if(!all(names(evidence)%in%names(bn))){
    unkownVariable = names(evidence)[which(!(names(evidence)%in%names(bn)))]
    stop("Some evidence names do not match any model names: ", paste(unkownVariable, collapse = ", "))
  }
  
  # Check values
  for(i in 1:length(bn)){
    if(bn[[i]]$node%in%names(evidence)){
      domain = bn[[i]]$functions[[1]]$Fx$Domain
      evi = evidence[[bn[[i]]$node]]
      if(between(evi, domain)==FALSE){
        stop("The evidence value of ", bn[[i]]$node," is outside its domain: ", paste0('[', domain[1], ', ', domain[2], ']'))
      }
    }
  }
  #******************************
  #*
  #****** nodo observado ******
  # obs = which(names(bn)%in%names(evidence))
  obs = names(bn)[names(bn)%in%names(evidence)]
  for(i in 1:length(obs)){
    evi = evidence[obs[i]]
    fx = bn[[obs[i]]]$functions
    # evaluar mop con valor observado
    for(k in 1:length(fx)){
      bn[[obs[i]]]$functions[[k]]$Fx$Function = eval.motbf(fx[[k]]$Fx, evi)
    }
  }
  
  #****** nodo hijo de nodos observados ****
  # variables involved in each MOP (node and parents). Scope of each factor
  factorScope = lapply(bn, attr, "scope")
  
  # nodos padres de la RB
  parents = lapply(factorScope, function(x){x[-length(x)]})
  
  # nodo cuyo padre es observado
  # (intervalos cerrados por la derecha: a<x<=b)
  child = which(detect_string(parents, names(evidence))==TRUE)
  
  if(length(child)>0){
    for(i in 1:length(child)){
      obs_parents = names(evidence)[names(evidence)%in%parents[[child[i]]]]
      fx = bn[[child[i]]]$functions
      
      # encontrar los mop validos para la observacion dada
      # keep = which(sapply(sapply(fx, '[',1), keep_fx, evi = as.data.frame(evidence)) == TRUE)
      parentIntervals = sapply(fx, '[',1)
      evi = as.data.frame(evidence)
      keep = which(keep_fx(evi, parentIntervals))
      setNode = bn[[child[i]]]$functions[keep]
      

      bn[[child[i]]]$functions = setNode
    }
  }
  
  print(bn)
  return(invisible(bn))
}

# for networks of class "motbf.ve"
# returns "mutilated network"
setEvidence_VE <- function(bn, evidence){
  
  # #--------- COMPROBAR EVIDENCIA -----------#
  # # Check names
  # if(!all(names(evidence)%in%names(bn))){
  #   unkownVariable = names(evidence)[which(!(names(evidence)%in%names(bn)))]
  #   stop("Some evidence names do not match any model names: ", paste(unkownVariable, collapse = ", "))
  # }
  # 
  # # Check values
  # 
  # for(i in 1:length(bn)){
  #   node = bn[[i]]$node
  #   if(node%in%names(evidence)){
  #     # domain = bn[[i]]$functions[[1]]$Fx$Domain
  #     domain = levels(bn)[[node]]
  #     evi = evidence[[node]]
  #     if(bn[[i]]$type == 'Discrete'){
  #       if(!any(evi == domain)){
  #         stop("The evidence value of ", node," is not included in its set of possible values: ", paste0('{', paste(domain, collapse = ", "), '}'))
  #       }
  #     }else{
  #       if(between(evi, domain)==FALSE){
  #         stop("The evidence value ", evi," of ", node," is outside its domain: ", paste0('[', domain[1], ', ', domain[2], ']'))
  #       }
  #     }
  #     
  #   }
  # }
  
  #------------------------------------------#
  #---------- nodo observado ----------------#
  #------------------------------------------#
  # obs = which(names(bn)%in%names(evidence))
  
  obs = names(bn)[names(bn)%in%names(evidence)]
  
  # variables involved in each MOP (node and parents). Scope of each factor
  factorScope = lapply(bn, attr, "scope")
  
  for(i in 1:length(obs)){
    observedNode = bn[[obs[i]]]
    evi = evidence[obs[i]]
    if(observedNode$type=='Discrete'){ # discrete observed node ----#
      # Fing all tables where the observe node appears (i.e., its own and its children's)
      mop_ind = which(detect_string(factorScope, observedNode$node)==TRUE)
      mop_ind
      
      
      for(k in 1:length(mop_ind)){
        # subset the table compatible with evidence
        keep = which(evi ==bn[[mop_ind[k]]]$functions[, observedNode$node, drop = FALSE])
        bn[[mop_ind[k]]]$functions = bn[[mop_ind[k]]]$functions[keep,]
      }
      # sapply(bn[mop_ind], '[[', 'functions')
    }else{ # continuous observed node ----#
      fx = observedNode$functions
      # evaluar mop con valor observado
      
      n = ncol(fx)
      for(k in 1:nrow(fx)){# each row of fx corresponds to a function
        # the function is stored in the last column of object 'fx'
        bn[[obs[i]]]$functions[[k,n]] = new_mop(eval.motbf(fx[k,n][[1]], evi))
      }
      
      # Find children of this observed node 
      parents = lapply(factorScope, function(x){x[-length(x)]})
      child = which(detect_string(parents, observedNode$node)==TRUE)
      
      if(length(child)>0){
        for(k in 1:length(child)){
          fx = bn[[child[k]]]$functions
          
          # intervals of observed parent in child's node
          A = fx[[observedNode$node]]
          
          # reformat for keep_fx()
          A = lapply(A, function(x){
            y = data.frame(x)
            colnames(y) = observedNode$node
            return(y)
          }) 
          
          # (intervalos cerrados por la derecha: a<x<=b)
          keep = which(keep_fx(as.data.frame(evi), A))
          bn[[child[k]]]$functions = fx[keep,]
        }
      }
    }
  }

  lev = levels(bn)
  evidence = evidence[match(names(lev[names(lev)%in%names(evidence)]), names(evidence))]
  lev[which(names(lev)%in%names(evidence))]  = evidence
  
  attr(bn, 'levels') = lev
  # print(sapply(bn, '[[', 'functions'))
  
  # Incompatible evidence under given model
  checkInconsistentEvidence(bn)
  
  return(invisible(bn))
}

checkInconsistentEvidence <- function(bn) {

  ## combination of evidence does not exist
  # this produces lists with empty data.frames (nrow = 0)
  A = lapply(bn, '[[', 'functions')
  
  if(any(sapply(A, nrow) == 0)){
    stop('inconsistent evidence has been attempted')
  }
  
  ## combination of evidence produces probability = 0
  # this occurs in discrete nodes
  
  B = sapply(A, function(x){x[,c(length(x)), drop = FALSE]})
  
  discrete.nodes = which(sapply(bn, '[[', 'type') == 'Discrete')
  if(any(sapply(B[discrete.nodes], function(x){all(x == 0)}))){
    stop('inconsistent evidence has been attempted')
  }
}
#**************************************
keep_fx <- function(evi, parentIntervals){

  res_fin <- c()
  
  minGlobal = min(sapply(parentIntervals, function(x){sapply(x, min)}))
  
  for(k in 1:length(parentIntervals)){
    obs = parentIntervals[[k]]
    n = colnames(evi)[which(colnames(evi)%in%colnames(obs))]
    
    res <- c()
    for(i in 1:length(n)){
      intervals = obs[,n[i]]
      
      x = evi[,n[i]]
      
      lower = intervals[1]
      upper = intervals[2]
      
      if(minGlobal == lower){
        lower = lower-0.1
      }
      ans = lower<x & x<=upper
      # if(k == 1){
      #   if(lower<=x & x<=upper){
      #     ans = TRUE
      #   }else{
      #     ans = FALSE
      #   }
      # }else{
      #   if(lower<x & x<=upper){
      #     ans = TRUE
      #   }else{
      #     ans = FALSE
      #   }
      # }
      
      res[i] <- ans
    }
    res_fin[k] <- all(res)
  }
  
  return(res_fin)
}


# between <- function(evi, intervals){
#   
#   lower = intervals[1]
#   upper = intervals[2]
#   
#   x = evi
#   if((lower <= x | abs(lower-x)<10^-16) & (x <= upper | abs(upper-x)<10^-16)){
#     ans = TRUE
#   }else{
#     ans = FALSE
#   }
#   return(ans)
# }







detect_string <- function(list, string){ 
  
  # lapply(list, function(x) { 
  #   if(TRUE %in% str_detect(x, paste(string, collapse = "|"))){ 
  #     TRUE 
  #   } else {
  #     FALSE
  #   } 
  # })
  
  lapply(list, function(x){
    string%in%x
  })
} 






normalizeMOP <- function(fx){

  f = fx$functions
  dom = f[,1]
  
  p = ncol(f)
  
  # detect if variable is discrete or continuous
  if(is.numeric(unlist(dom))){# continuous variable
    
    f = f[,p]
    x = c()
    
    # Definite integral in each split
    for(i in 1:length(f)){
      x[i] = integrate.motbf(f[[i]], f[[i]]$Domain[1,getMotbfVar(f[[i]])], f[[i]]$Domain[2,getMotbfVar(f[[i]])])
    }
    
    # normalization constant
    k = 1/sum(x)
    g = list()
    for(i in 1:length(f)){
      g[[i]] = multiplyMOPbyConstant(f[[i]], k)
    }
    
    
    # check if g(x) integrates to 1
    y = c()
    for(i in 1:length(g)){
      y[i] = integrate.motbf(g[[i]], g[[i]]$Domain[1,getMotbfVar(g[[i]])], g[[i]]$Domain[2, getMotbfVar(g[[i]])])
    }
    if(abs(1-sum(y))<10^-6){ # error smaller than 1e-6
      return(g)
    }else{
      stop('The MOP density does not integrate to 1')
    }
  }else{ # discrete variable
    pmf = f[[p]]
    # x = sapply(lapply(pmf, '[[', 'Function'), as.numeric)
    x = sapply(lapply(pmf, '[[', 1), as.numeric)
    k = 1/sum(x)
    
    g = list()
    for(i in 1:length(pmf)){
      g[[i]] = multiplyMOPbyConstant(pmf[[i]], k)
    }
    
    # check if g(x) sums up to 1
    y = sapply(lapply(g, '[[', 1), as.numeric)
    
    if(abs(1-sum(y))<10^-6){ # error smaller than 1e-6
      f[[p]] = g
      colnames(f)[2] = 'Marginal_distribution'
      return(f)
    }else{
      stop('The PMF does not sum up to 1')
    }
  }
  

}


# .-ESTA FUNCION NO SE USA ----
# divideMOPbyConstant = function(poly, cte, range = NULL){
#   
#   
#   st = splitMOP(poly)
#   gr = st$Exponents
#   cf = st$Coefficients/cte
#   
#   s = ifelse(sign(cf)>=0, "+", "")
#   s[1] = ''
#   t4 = paste0(s,cf, gr, collapse = "")
#   
#   if(is.null(range)){
#     range = poly$Domain
#   }
#   if(is.motbf(poly)){
#     t5 <- list(Function=noquote(t4), Domain= range, Subclass="mop")
#     output <- motbf(t5)
#   }else if(is.jointmotbf(poly)){
#     t5 <- list(Function=noquote(t4), Domain= range)
#     output <- jointmotbf(t5)
#   }
#   
#   
#   return (output)
# }


multiplyMOPbyConstant = function(poly, cte, range = NULL){
  
  
  st = splitMOP(poly)
  gr = st$Exponents
  cf = st$Coefficients*cte
  
  s = ifelse(sign(cf)>=0, "+", "")
  s[1] = ''
  t4 = paste0(s,cf, gr, collapse = "")
  
  if(is.null(range)){
    range = tryCatch({poly$Domain}, error=function(e){NULL})
    # range = poly$Domain
  }
  if(is.univmotbf(poly)){
    t5 <- list(Function=noquote(t4), Domain= range, Subclass="mop")
    output <- new_mop(t5)
  }else if(is.jointmotbf(poly)){
    t5 <- list(Function=noquote(t4), Domain= range)
    output <- new_jointmotbf(t5)
  }
  
  
  return (output)
}

getTopoOrder <- function(x){
  if(is.motbf_fit(x)){
    x = lapply(x, function(x){
      pa = x$parents
      node = x$node
      c(pa, node)
    })
  }
  vars = unique(unlist(x))
  m = matrix(rep(0,length(vars)^2), ncol = length(vars), dimnames = list(vars, vars))
  m
  for(i in 1:length(x)){
    node = x[[i]]
    child = node[length(node)]
    parent = node[-length(node)]
    
    r = which(rownames(m) %in% parent)
    c = which(colnames(m) %in% child)
    
    m[r,c] = 1
  }
  
  to = colnames(m)[topOrder(m)]
  return(to)
}


#' Export discrete motbf to bnlearn format
#' 
#' Export a discrete BN created with \link{motbf.fit} to an object of class \code{"bn"} of package bnlearn.
#' 
#' @param bn An object of class \code{motbf_fit}, obtained from function \link{motbf.fit}. Only discrete networks are valid.
#' @return An object of class \code{bn}.
#' @importFrom bnlearn custom.fit
#' @export
motbf2bnlearn = function(bn){
  if(!all(sapply(bn, '[[','type') == 'Discrete')){
    stop('The bn object must be discrete')
  }
  # if(!("motbf.ve" %in%class(bn))){
  #   bn = getFormatedBN_VE(bn)
  # }
  
  dag = getDAG(bn)
  
  variables = names(bn)
  
  saveNodes = list()

  for(i in 1:length(variables)){
    X = variables[i]
    # if(X == 'SATURACION'){
    #   browser()
    # }

    # cptX_dimnames = levels(bn)[c(X, bnlearn::parents(dag, X))]
    n = ncol(bn[[X]]$functions)
    p_dist = bn[[X]]$functions[,-c(n,(n-1)), drop = FALSE]
    
    # p = bnlearn::parents(dag, X)
    p = colnames(p_dist)
    if(length(p)>1){
      p = p[length(p):1]
    }
    cptX_dimnames = levels(bn)[c(X, p)]
    
    cptX_dim = sapply(cptX_dimnames, length)
    
    cptX = bn[[X]]$functions[,ncol(bn[[X]]$functions)]
    
    if(nrow(p_dist)!=prod(cptX_dim)){
      combinations = expand.grid(levels(bn)[names(cptX_dimnames)])[,-1, drop = FALSE]
      if(ncol(combinations)>1){
        combinations = combinations[,ncol(combinations):1]
      }
      
      A = apply(p_dist, 1, paste, collapse = ',')
      B = apply(combinations, 1, paste, collapse = ',')
      
      # suppressWarnings({
      #   cptX[which(A!=B)] = 1/length(which(A!=B))
      # })
      
      # which combinations are not observed in the data
      combNotObs = unique(B[which(!(B%in%A))])
      if(length(combNotObs)>0){
        # get uniform distribution for combinations not observed in the data
        unifDistNotObs = 1/length(which(B==combNotObs[1])) 
        # logic vector: T if combinations are observed; F if not
        cptLogic = B%in%A
        # keep CPT where combinations are observed
        cptLogic[which(cptLogic == TRUE)] = cptX
        # add uniform distribution where combinations are not observed
        cptLogic[!(B%in%A)] = unifDistNotObs
        
        cptX = cptLogic
      }
    }
    
    dim(cptX) = cptX_dim
    dimnames(cptX) = cptX_dimnames
    saveNodes[[i]] = cptX
  }
  
  names(saveNodes) = variables
  dfit = custom.fit(dag, dist = saveNodes)
  
  return(dfit)
}
#----


#' Export discrete motbf to grain format
#' 
#' Export a discrete BN created with \link{motbf.fit} to an object of class \code{"grain"} of package gRain
#' 
#' @param bn An object of class \code{motbf_fit}, obtained from function \link{motbf.fit}. Only discrete networks are valid.
#' @return An object of class \code{grain}.
#' @importFrom bnlearn as.grain
#' @export
motbf2grain = function(bn){
  if(!all(sapply(bn, '[[','type') == 'Discrete')){
    stop('The bn object must be discrete')
  }
  
  asbnfit = motbf2bnlearn(bn)
  asgrain = as.grain(asbnfit)

  return(asgrain)
}