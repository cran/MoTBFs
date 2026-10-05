#' Fitting MoTBFs TAN models
#'
#' Perform a TAN model of class MoTBF based on maximizing the Mutual Information.
#' 
#' 
#' @name MOPTAN
#' @rdname MOPTAN
#' @param data A \code{data.frame} containing the variables, which can contain 
#' continuous and discrete variables.
#' @param target A \code{character} string indicating the name of the target or 
#' class variable. Target must be a column of \code{data} and it can be a 
#' continuous or discrete variable.
#' @param fit.args A \code{list} a list containing optional arguments used to 
#' fit the models. 
#' These arguments must be those accepted by function \link{motbf.fit},
#' i.e., 'numIntervals' (4), 'POTENTIAL_TYPE' ('MOP'), 'maxParam' (7),
#' 's' (NULL), 'priorData' (NULL) or 'scale' (TRUE).
#' If fit.args is left NULL, the default values (in brackets) for those 
#' arguments will be used.
#' @param root A \code{character} string indicating the label of the root
#' predictor variable of TAN model.
#' @param all A \code{logical} flag. If \code{TRUE}, the function return a 
#' \code{list} which contains two elements: the TAN model and the mutual 
#' information used to compute the maximum spanning tree. 
#' Defaults to \code{FALSE}.
#' @param mutualInfoCond A numeric matrix indicating the estimation of the 
#' mutual information coefficientes for Chow-Liu- algorithm in TAN. If it is 
#' \code{NULL}, it is computed using \link{mutual_information_tan} function.
#' @param parallel A \code{logical} flag. If \code{TRUE}, computation runs in parallel 
#'   using \code{foreach} and \code{doParallel}. Defaults to \code{FALSE}.
#' @importFrom parallel detectCores mclapply makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' 
#' @details
#' The main function, \code{fit_tan()}, fits a MoTBF Tree Augmented Naive Bayes
#' model using the specified \code{data}.
#' 
#' 
#' @return 
#' The main function, \code{fit_tan()}, returns an object of class \code{"motbf_fit"}.
#' When \code{all=TRUE}, it returns a list with the Bayesian network and the mutual
#' information matrix used to compute the maximun spanning tree. 
#' 
#' Function \code{mutual_information_tan()} returns a symmetric numeric 
#' \code{matrix} of dimensions \eqn{k \times k}{k x k}, where \eqn{k} is the number of 
#' predictor variables. Row and column names correspond to the predictor 
#' variables, and the entries contain the estimated conditional mutual 
#' information values. This matrix is used to compute the TAN model.
#' 
#' @examples 
#' \donttest{
#' data = iris
#' data$Species = as.factor(data$Species)
#' # Fit TAN model for classification
#' tan = fit_tan("Species",data)
#' }
#' 
## Funcion para aprender un tan-------------------------------------
#' @export
fit_tan <- function(target, data,fit.args=NULL, root = NULL, all = FALSE,
                    mutualInfoCond = NULL, parallel = FALSE){
  
  if(is.null(target)){
    stop(paste0("Variable ", target, " must be indicated")) 
  }
  
  if(!(target%in%colnames(data))){
    stop(paste0("Variable ", target, " must be a name of a variable of data")) 
  }
  if(is.null(fit.args)){
    fit.args=fit.args.null(fit.args)
  }
  
  data = check_data(data)
  if(!is.null(fit.args$priorData)){
    priorData <- check_data(fit.args$priorData)
    if(scale){
      priorData_original <- fit.args$priorData
      fit.args$priorData <- rescale_data(fit.args$priorData)
    }
  }
  scale = fit.args$scale
  if(scale){
    data_original = data
    data = rescale_data(data)
  }
  fit.args$scale=FALSE
  # browser()
  # Obtencion de la informacion mutua y distribuciones-------------------
  # Obtencion de la estructura--------------------------------------
  
  if(is.null(mutualInfoCond)){
    mutualInfoCond = mutual_information_tan(data,target = target,fit.args = fit.args,
                                            parallel = parallel)
    
  }
  if(is.null(root)){
    root = fitroot(mutualInfoCond)
  }
  # Calculo del grafo
  # Matriz de adjacencia
  adj_mat <- prim_maximal(mutualInfoCond,root = root)
  # nombre de las variables
  varTAN <- c(row.names(adj_mat),target)
  # anadimos la clase
  adj_mat <- rbind(adj_mat,1)# Aristas del target a todos los predictores
  adj_mat <- cbind(adj_mat,0)
  rownames(adj_mat)<- varTAN
  colnames(adj_mat)<- varTAN
  
  
  varPred = setdiff(colnames(data),target)
  # var_i = colnames(data)[1]
  
  parents = sapply(varPred,function(var_i){
    names(which(adj_mat[,var_i]==1))
  })
  if(parallel){
    
    cores = detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    # browser()
    var_i <- NULL
    distributions = foreach(var_i = varPred,.packages = "MoTBFs",
                            .inorder = FALSE)%dopar%{
                              distX_cond(data,var_i,parents[[var_i]],fit.args)
                            }
    stopCluster(cl)
  }else{
    # browser()
    distributions = lapply(varPred,function(var_i){
      distX_cond(data,var_i,parents[[var_i]],fit.args)
    })
  }
  # distribucion de la clase
  distTarget = distVar(data = data,target,fit.args)
  distributions = append(distributions,list(distTarget))
  
  bn = getFormatedBN(distributions)
  # dag = getDAG(bn)
  # graphviz.plot(dag)
  if(scale){
    bn = rescale_motbf_fit(bn, POTENTIAL_TYPE = fit.args$POTENTIAL_TYPE, data = data_original)
  }
  if(all){# Devolver informacion mutua y distribicones tambien
    return(list(bn=bn,mutualInfoCond=mutualInfoCond))
  }else{
    return(bn)
  }
}
## Funcion para calcular la matriz de IM dada la clase-------------------
#' @importFrom parallel detectCores mclapply makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @export
#' @rdname MOPTAN
mutual_information_tan = function(data,target,fit.args = NULL,parallel = FALSE){
  # browser()
  # source("mutualInformation.R")
  data = check_data(data)
  fit.args=fit.args.null(fit.args)
  # Seleccionamos las variables
  var_i <- i <- NULL
  varNames = colnames(data)
  # variables predictoras
  varPred = varNames[varNames!=target]
  
  
  # Inicializamos la matriz para MI
  # Inicializamos la matriz para MI
  MI = matrix(0,nrow = length(varPred),ncol = length(varPred),
              dimnames = list(varPred,varPred))
  ff = function(i,X,target,data,distC,distX_C,distributions,fit.args){
    Y = varPred[i]
    # browser()
    if(X==Y){
      return(list(MI=0,name = Y))
    }else{
      distY_C = distributions[[paste0(Y,"|",target)]]
      # print(distY_C)
      # exists("cond_mut_information_cont")
      result = tryCatch(cond_mut_information_compute(data = data,className = target,
                                                     # varNames = c(X,Y,target),
                                                     fit.args = fit.args,
                                                     distC = distC,distY_C = distY_C,
                                                     distX_C = distX_C),
                        error = function(e) list(MI=0))
      return(list(MI = result$MI,name = Y))
    }
  }
  if(parallel){
    time = Sys.time()

    # setup parallel backend to use many processors
    cores = detectCores()
    # cl <- makeCluster(cores[1]-1) #not to overload your computer
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    result = foreach(var_i = varPred,.packages = "MoTBFs")%dopar%{
      distX_C = distX_cond(data = data,var_i,target,
                           fit.args = fit.args)
      distX_C = getFormatedBN(list(distX_C))[[1]]$functions
      colnames(distX_C)[3]=paste0("CPD_",var_i,"|",target)
      list(distX_C,paste0(var_i,"|",target))
    }
    # Aprendizaje distribuciones X|C
    distributions = lapply(result,"[[",1)
    names(distributions)=lapply(result,"[[",2)
    distC = distVar(data[,target,drop = FALSE],target,fit.args)
    ## Nos quedamos solo con la CPD
    distC =  getFormatedBN(list(distC))[[1]]$functions
    # Calcular la MI
    foreach(X = varPred,
            .packages = c("MoTBFs"))%do%{
              # variable que es condicionada
              # Aprendizaje de la distribucion var_i|var_j,C
              distX_C = distributions[[paste0(X,"|",target)]]
              
              
              result = foreach(i = 1:length(varPred),
                               .packages = c("MoTBFs"))%dopar%{
                                 ff(i,X,target,data,distC,distX_C,
                                    distributions,fit.args)
                               }
              
              
              MIs = sapply(result, "[[", "MI")
              namesDist = sapply(result, "[[", "name")
              # Padres por filas e hijos por columnas
              MI[sapply(namesDist, "[", 1),X]=MIs
            }
    stopCluster(cl)
    Sys.time()-time
  }else{
    # Distribucion de la variable clase y de X_i|C--------------------
    distributions = lapply(varPred, function(X){
      distX_C = distX_cond(data = data,varX = X,
                           varsCond = target,fit.args = fit.args)
      distX_C = getFormatedBN(list(distX_C))[[1]]$functions
      colnames(distX_C)[3]=paste0("CPD_",X,"|",target)
      return(distX_C)
    })
    names(distributions)=paste0(varPred,"|",target)
    distC = distVar(data,varName = target,fit.args = fit.args)
    ## Nos quedamos solo con la CPD
    distC =  getFormatedBN(list(distC))[[1]]$functions
    # X = varPred[3]
    # nDist = length(distributions)
    # browser()
    for(X in varPred){
      # variable que es condicionada
      # Aprendizaje de la distribucion var_i|var_j,C
      distX_C = distributions[[paste0(X,"|",target)]]
      # Y = varPred[2]
      for(Y in varPred){
        # iter = iter+1
        if(X==Y){
          next
        }
        distY_C = distributions[[paste0(Y,"|",target)]]
        result=cond_mut_information_compute(data = data[,c(X,Y,target)],
                                            className = target,
                                            fit.args = fit.args,
                                            distC = distC,distY_C = distY_C,
                                            distX_C = distX_C)
        MI[Y,X]=result$MI
      }
    }
  }
  # Hacemos la media para considerar un grafo no dirigido
  MI = (MI +t(MI))/2
  return(MI)
}



# Funciones de control------------------------------------------

#' @noRd
fit.args.null = function(fit.args){
  if(is.null(fit.args$numIntervals)){
    fit.args$numIntervals=4
  }
  if(is.null(fit.args$POTENTIAL_TYPE)){
    fit.args$POTENTIAL_TYPE="MOP"
  }
  if(is.null(fit.args$scale)){
    fit.args$scale=FALSE
  }
  if(is.null(fit.args$maxParam)){
    fit.args$maxParam = 7
  }
  return(fit.args)
}

# Funciones para calcular distribuciones---------------------------------------
#' @noRd
distX_cond = function(data,varX, varsCond, fit.args =NULL){
  # browser()
  varType = ifelse(is.numeric(data[,varX]),"Continuous","Discrete")
  distX_C = do.call(conditionalMethod,
                    append(fit.args,
                           list(data=data,nameParents=varsCond,
                                nameChild=varX)))
  distX_C = list(Child=varX,functions=distX_C,varType = varType)
  
  return(distX_C)
}


#' @noRd
distVar = function(data,varName,fit.args){
  fit.args$s = NULL
  fit.args$priorData = NULL
  if(is.numeric(data[,varName])){
    fit.args$x=data[,varName]
    fit.args$numIntervals=NULL
    distVar = do.call(univMoTBF,fit.args)
    distVar = list(Child=varName,functions=list(Px = distVar),# Para que este en el formato correcto
                   varType = "Continuous")
  }else{
    distVar = probDiscreteVariable(data[,varName])
    # Almacenamos la distribucion para devolverla
    distVar = list(Child=varName,functions=list(distVar),varType = "Discrete")
  }
  return(distVar)
}
# Funciones para integrar MoTBFs-----------------------------------------------
# Funcion para integrar en X
#' @noRd
integrateX = function(CPD_X_,CPD_X,X){
  # browser()
  integrateX = c()
  for(i in 1:length(X)){
    domain = unlist(X[i])
    # tranformamos en funcion
    CPD_X_i = as.function(CPD_X_[i][[1]])
    CPD_Xi = as.function(CPD_X[i][[1]])
    # Calculo del integrando
    # integrando = function(x){
    #   xx = seq(domain[1],domain[2],10^-4)
    #   minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
    #   return(CPD_X_i(x)*(log(CPD_X_i(x)/minn)-log(CPD_Xi(x)/minn)))
    #   # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    # }
    # integrateX[i]=integrate(integrando,domain[1],domain[2])$value
    integrando1 = function(x){
      # xx = seq(domain[1],domain[2],10^-4)
      # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
      return(CPD_X_i(x)*log(CPD_X_i(x)))
      # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    }
    integrando2 = function(x){
      # xx = seq(domain[1],domain[2],10^-4)
      # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
      return(CPD_X_i(x)*log(CPD_Xi(x)))
      # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
    }
    integrateX[i]=tryCatch(integrate(integrando1,domain[1],domain[2])$value-
                             integrate(integrando2,domain[1],domain[2])$value,
                           error = function(e) {
                             cat("Error:", e$message, "\n")
                             0},warning = function(w) {
                               cat("Advertencia:", w$message, "\n")
                               0})
  }
  return(integrateX)
}

#' @noRd
integrateVar = function(CPD,domain){
  
  integrateVar=c()
  i = 1
  for(i in 1:length(domain)){
    CPD_i = CPD[i][[1]]
    domain_i = unlist(domain[i])
    integrateVar[i] = integrate.motbf(CPD_i,domain_i[1],domain_i[2])
  }
  return(integrateVar)
}

# Informacion mutua para dos variables discretas-----------------

# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introduccir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X la variable continua e Y la variable discreta. Todas aprendidas con 
# conditionalMethod


## Funcion para calcular el valor (auxiliar no se exporta)----------------
#' @noRd
mut_information_compute = function(data, nameVars=NULL, fit.args =NULL, all = FALSE,
                                   distY = NULL,distX_Y=NULL,distX = NULL){

  # Control nameVars
  if(!is.null(nameVars)){
    data = data[,nameVars]
    if(is.null(data)){
      stop(paste0("data must have",nameVars,"as colnames",collapse = " "))
    }
  }
  if(length(data)!=2){
    stop("data must have 2 variables")
  }else{
    nameVars = colnames(data)
  }
  # browser()
  X = nameVars[1]
  Y = nameVars[2]
  distributions = list()
  nDist=0
  ## Distribucion de X|Y------------------------------------
  if(is.null(distX_Y)){
    # distX_Y = conditionalMethod(data = data, nameParents = c("Y"),
    #                                 nameChild = "X",numIntervals = 4,
    #                                 POTENTIAL_TYPE = "MOP",scale = scale)
    distX_Y = distX_cond(data = data,varX = X,
                         varsCond = Y,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist = nDist+1
    distributions[[nDist]] = distX_Y
    names(distributions)[nDist]=paste0(X,"|",Y)
    ## Nos quedamos solo con la CPD
    distX_Y = getFormatedBN(list(distX_Y))[[1]]$functions
  }else{
    ## Nos quedamos solo con la CPD
    distX_Y = getFormatedBN(distX_Y)[[1]]$functions
  }
  colnames(distX_Y)[3]=paste0("CPD_",X,"|",Y)
  
  # Distribucion de Y----------------------------
  if(is.null(distY)){
    distY = distVar(data,varName = Y,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distY
    names(distributions)[nDist]=Y
    ## Nos quedamos solo con la CPD
    distY =  getFormatedBN(list(distY))[[1]]$functions
  }else{
    distY =  getFormatedBN(distY)[[1]]$functions
  }
  
  # Distribucion de X----------------------------
  if(is.null(distX)){
    distX = distVar(data,varName = X,fit.args = fit.args)
    # Almacenamos la distribucion para devolverla
    nDist=nDist+1
    distributions[[nDist]] = distX
    names(distributions)[nDist]=X
    ## Nos quedamos solo con la CPD
    distX =  getFormatedBN(list(distX))[[1]]$functions
  }else{
    distX =  getFormatedBN(distX)[[1]]$functions
  }
  
  
  # Calculo de la MI
  numLevels = sapply(data,nlevels)
  # browser()
  if(all(numLevels==0)){
    # Variables continuas
    # Bucle que pasa por cada particion
    MI = c()
    i = 1
    px = as.function(distX[[2]][[1]])
    py = as.function(distY[[2]][[1]])
    for(i in 1:nrow(distX_Y)){
      # Seleccionamos la distribucion X|Y \in y_i
      px_yi = as.function(distX_Y[[i,paste0("CPD_",X,"|",Y)]])
      # DOminio de la variable X|Y\in I_i
      domain = sort(distX_Y[[i,X]])
      # Intervalo de la particion de la variable Y
      intervals = sort(distX_Y[[i,Y]])
      
      # Calculo de la primera integral (X|Y \in y_i)
      integrando = function(x){
        return(px_yi(x)*log(px_yi(x)/px(x)))
      }
      integral = integrate(integrando,domain[1],domain[2])$value
      # Calculo integral Y*Integral1 si Y \in y_i
      MI = c(MI,integrate(py,intervals[1],intervals[2])$value*integral)
      MI = sum(MI)
    }
  }else if(all(numLevels>1)){
    # Variables discretas
    # Arbol para los calculos
    arbol = merge.data.frame(distX_Y,distX)
    arbol = merge.data.frame(arbol,distY)
    
    # Calculo de la informacion para cada rama
    arbol$sumX = arbol[,paste0("CPD_",X,"|",Y)]*
      (log(arbol[,paste0("CPD_",X,"|",Y)])-log(arbol[,paste0("CPD_",X)]))
    # Compute I(X,Y)
    MI = sum(arbol$sumX*arbol[,paste0("CPD_",Y)])
  }else{
    # Variable continua y discreta
    # Condiciona la variable discreta:
    if(!is.numeric(data[,Y])){
      # Construimos el arbol para los calculos
      arbol = cbind(distX_Y,distX[,2,drop=FALSE])
      # Incluimos la probabilidad de Y
      arbol = merge.data.frame(arbol,distY)
      # Calculamos las diferentes integrales en X
      arbol$integrateX = integrateX(arbol[,paste0("CPD_",X,"|",Y)],
                                    arbol[,paste0("CPD_",X)],
                                    arbol[,X])
      MI = sum(arbol[,paste0("CPD_",Y)]*arbol$integrateX)
    }else{
      # Condiciona la variable continua
      # Construimos el arbol para los calculos
      arbol = merge.data.frame(distX_Y,distX)
      # Incluimos la probabilidad de Y
      arbol = merge.data.frame(arbol,distY[,2,drop=FALSE])
      # Calculamos las diferentes integrales en X
      arbol$sumX = log(arbol[[paste0("CPD_",X,"|",Y)]]/arbol[[paste0("CPD_",X)]])*
        arbol[[paste0("CPD_",X,"|",Y)]]
      arbol$integrateY = integrateVar(arbol[[paste0("CPD_",Y)]],domain = arbol[[Y]])
      MI = sum(arbol$sumX*arbol$integrateY)
    }
  }
  if(all){
    return(list(MI = MI, distributions = distributions))
  }
  return(MI)
}

## Funcion para calcular la MI entre dos variables (si se exporta)-------
#' @noRd
mut_information = function(data, nameVars=NULL, fit.args =NULL, all = FALSE){
  # Chequeamos los argumentos para fijar las Mops------------
   data = check_data(data)
  fit.args = fit.args.null(fit.args)
  # Control nameVars
  if(!is.null(nameVars)){
    data = data[,nameVars]
    if(is.null(data)){
      stop(paste0("data must have",nameVars,"as colnames",collapse = " "))
    }
  }
  if(length(data)!=2){
    stop("data must have 2 variables or nameVars must have two values")
  }else{
    nameVars = colnames(data)
  }
  # Condiciona Y
  MI1 = mut_information_compute(data,nameVars = nameVars,fit.args,all = TRUE)
  # Condiciona X
  MI2 = mut_information_compute(data,nameVars = nameVars[2:1],
                                fit.args = fit.args,all=FALSE,
                                distX = MI1$distributions[2],
                                distY = MI1$distributions[3])
  return((MI1$MI+MI2)/2)
}


 
# # Variables continuas y discretas
# 
# 
# Y = as.factor(c(rep("y1",40),rep("y2",60)))
# set.seed(1252)
# X = rbeta(40,2,4)
# set.seed(1262)
# X = c(X,rbeta(60,5,1))
# data = data.frame(X,Y)
# mut_information(data)
# 
# # Variables continuas
# set.seed(1252)
# X = rbeta(100,4,4)
# set.seed(1252)
# Y = rbeta(100,3,3)
# data = data.frame(X,Y)
# mut_information(data)

# Informacion mutua entre variables dada otra variable------



## Funcion para no exportar------------------------------

# fit.args =  Argumentos para la calcular las distribuciones MoTBFs
# Vease argumentos de ConditionalMethod. Si es null, emplea los
# argumentos por defecto de esa funcion con ademas numIntervals=4,
# POTENTIAL_TYPE = "MOP",scale = F
# Se puede introduccir la distribucion de C obtenida empleando la
# la funcion probDiscreteVariable()
# se puede introduccir tambien las distribuciones de Y|C, de X|C y de X| Y,C,
# con X la variable continua e Y la variable discreta. Todas aprendidas con 
# conditionalMethod

#' @noRd
cond_mut_information_compute = function(data,className, fit.args,
                                        distC = NULL,distY_C=NULL,
                                        distX_YC=NULL,distX_C=NULL){
  
  # Eliminamos la variable clase del data.frame
  data2 = data[,colnames(data)!=className]
  
  # browser()
  # Calculamos las distribuciones
  X = colnames(data2)[1]
  X_discrete = !is.numeric(data[,X])
  Y = colnames(data2)[2]# Y es la variable que condiciona
  Y_discrete = !is.numeric(data[,Y])
  C_discrete = !is.numeric(data[,className])
  rm(data2)
  
  distributions = list()
  # Cuenta el numero nuevo de distribuciones aprendidas
  nDist = 0
  
  # Distribucion de X|C,Y---------------------------------
  # distX_YC = conditionalMethod(data,c("C","Y"),"X",4,"MOP",scale=F)
  if(is.null(distX_YC)){
    distX_YC =  distX_cond(data = data,varX = X,varsCond = c(Y,className),
                           fit.args = fit.args)
    ## Nos quedamos solo con la CPD
    distX_YC = getFormatedBN(list(distX_YC))[[1]]$functions
    colnames(distX_YC)[4]=paste0("CPD_",X,"|",Y,":",className)
  }
  # Almacenamos la distribucion para devolverla
  nDist=nDist+1
  distributions[[nDist]] = distX_YC
  parents=paste0(sort(c(className,Y)),collapse = ":")
  names(distributions)[nDist]=paste0(X,"|",Y,":",className)
  
  # Distribucion de X|C-----------------------------
  if(is.null(distX_C)){
    distX_C = distX_cond(data = data,varX = X,varsCond = c(className),
                         fit.args = fit.args)
    ## Nos quedamos solo con la CPD
    distX_C = getFormatedBN(list(distX_C))[[1]]$functions
    colnames(distX_C)[3]=paste0("CPD_",X,"|",className)
  }
  # Almacenamos la distribucion para devolverla
  nDist = nDist+1
  distributions[[nDist]] = distX_C
  names(distributions)[nDist]=paste0(X,"|",className)

  # Distribucion de C----------------------------
  if(is.null(distC)){
    distC = distVar(data,varName = className,fit.args = fit.args)
    ## Nos quedamos solo con la CPD
    distC =  getFormatedBN(list(distC))[[1]]$functions
  }
  # Almacenamos la distribucion para devolverla
  nDist = nDist+1
  distributions[[nDist]] = distC
  names(distributions)[nDist]=paste0(className)
  # Distribucion de Y|C--------------------------
  if(is.null(distY_C)){
    distY_C = distX_cond(data = data,varX = Y,varsCond = c(className),
                         fit.args = fit.args)
    ## Nos quedamos solo con la CPD
    distY_C = getFormatedBN(list(distY_C))[[1]]$functions
    colnames(distY_C)[3]=paste0("CPD_",Y,"|",className)
    
  }
  # Almacenamos la distribucion para devolverla
  nDist = nDist+1
  distributions[[nDist]] = distY_C
  names(distributions)[nDist]=paste0(Y,"|",className)
  # browser()
  # Calculo de la Informacion Mutua-------------------
  if(C_discrete){# variable clase discreta
    if(Y_discrete){# condiciona una discreta
      if(X_discrete){# discreta|discreta y clase discreta
        # Construimos el arbol para los calculos
        arbol = merge.data.frame(distX_YC,distX_C)
        # Incluimos la probabilidad de C
        arbol = merge.data.frame(arbol,distC)
        
        # Añadimos las probabilidades Y|C al arbol
        arbol = merge.data.frame(arbol,distY_C)
        
        # Calculo del log(p(X|YC)/p(X|C))
        arbol$log = arbol[,paste0("CPD_",X,"|",Y,":",className)]*
          (log(arbol[,paste0("CPD_",X,"|",Y,":",className)])-
             log(arbol[,paste0("CPD_",X,"|",className)]))
        
        # Calculamos la informacion mutua
        MI = sum(arbol$log*
                   arbol[,paste0("CPD_",Y,"|",className)]*
                   arbol[,paste0("CPD_",className)])
      }else{# continua|discreta y clase discreta
        # Construimos el arbol para los calculos
        arbol = merge.data.frame(distX_YC,distX_C)
        
        # Calculamos las diferentes integrales en X
        arbol$integrateX = integrateX(arbol[,paste0("CPD_",X,"|",Y,
                                                    ":",className)],
                                      arbol[,paste0("CPD_",X,"|",className)],
                                      arbol[,X])
        
        # Incluimos la probabilidad de C
        arbol = merge.data.frame(arbol,distC)
        
        # Añadimos las probabilidades Y|C al arbol
        arbol = merge.data.frame(arbol,distY_C)
        
        
        MI = sum(arbol$integrateX*
                   arbol[,paste0("CPD_",Y,"|",className)]*
                   arbol[,paste0("CPD_",className)])
      }
    }else{# condiciona una continua
      if(X_discrete){# Discreta|continua y clase discreta
        # Construimos el arbol para los calculos
        arbol = merge.data.frame(distX_YC,distX_C)
        # Incluimos la probabilidad de C
        arbol = merge.data.frame(arbol,distC)
        #Incluimos la distribucion de Y|C
        # no modificar el dominio de X
        arbol = merge.data.frame(arbol,distY_C[,-2])
        
        # p(x|yc)*(log(x|yc)-log(x|c))
        arbol$sumX = arbol[,paste0("CPD_",X,"|",Y,":",className)]*
          (log(arbol[,paste0("CPD_",X,"|",Y,":",className)])-
             log(arbol[,paste0("CPD_",X,"|",className)]))
        
        # Calculamos las diferentes integrales en Y
        arbol$integrateY = integrateVar(arbol[,paste0("CPD_",Y,"|",className)],
                                        arbol[,Y])
        # Calculo de la Informacion Mutua
        MI = sum(arbol$integrateY*arbol$sumX*arbol[,paste0("CPD_",className)])
      }else{# continua| continua y clase discreta
        # Construimos el arbol para los calculos
        arbol = merge.data.frame(distX_YC,distX_C)
        
        # Calculamos las diferentes integrales en X
        arbol$integrateX = integrateX(arbol[,paste0("CPD_",X,"|",Y,
                                                    ":",className)],
                                      arbol[,paste0("CPD_",X,"|",className)],
                                      arbol[,X])
        
        # Incluimos la probabilidad de C
        arbol = merge.data.frame(arbol,distC)
        
        # Añadimos las probabilidades Y|C al arbol
        arbol = merge.data.frame(arbol,distY_C[,-2])# No incluir el dominio de Y
        
        arbol$integrateY = integrateVar(arbol[,paste0("CPD_",Y,"|",className)],
                                        arbol[,Y])
        
        # Calculamos la informacion mutua
        MI = sum(arbol$integrateX*arbol$integrateY*arbol[,paste0("CPD_",className)])
      }
    }
  }else{# Variable clase continua
    MI = 0
    # Construimos el arbol para los calculos
    
    n_intervals = sapply(distributions,nrow)
    # posibles combinaciones de los intervalos de cada nodo
    intervals = lapply(n_intervals, function(x){seq(1:x)})
    combination.intervals = expand.grid(intervals)
    if(Y_discrete){# condiciona una discreta
      if(X_discrete){ # discreta|discreta y clase continua
        selected = list()
        k = 1
        i = 1
        j = 1
        for(i in 1:nrow(combination.intervals)){
          comb = combination.intervals[i,, drop = FALSE]
          
          
          for(j in 1:length(comb)){
            node = which(names(distributions) == names(comb[j])) # nodo
            selected[[j]] = distributions[[node]][unlist(comb[j]),] # intervalo dentro del nodo
            names(selected)[j] = names(comb)[j]
          }
          # sel2 = unlist(selected, recursive = F)
          
          # combined factor's domain
          domain = getFactorDomain(selected)
          
          if(!is.null(domain)){
            # sel2$DomainFactorParents = domain
            # res[[k]] = sel2
            
            # Obtenemos la funciones de probabilidad
            arbol = lapply(selected, function(x){x[[length(x)]]})
            names(arbol) = names(selected)
            
            arbol$X = arbol[[paste0(X,"|",Y,":",className)]]*log(arbol[[paste0(X,"|",Y,":",className)]])-
              arbol[[paste0(X,"|",Y,":",className)]]*log(arbol[[paste0(X,"|",className)]])
            arbol$integrateC = integrate.motbf(arbol[[paste0(className)]][[1]],
                                             domain[1,className],domain[2,className])
            
            MI = MI+arbol$integrateC*arbol$X*arbol[[paste0(Y,"|",className)]]
            
          }
          
        }
      }else{# continua|discreta y clase continua
        # comprobar si cada combinacion de intervalos se puede multiplicar
        # y obtener el dominio sobre el cual la funcion es valida
        
        selected = list()
        k = 1
        i = 1
        j = 1
        for(i in 1:nrow(combination.intervals)){
          comb = combination.intervals[i,, drop = FALSE]
          
          
          for(j in 1:length(comb)){
            node = which(names(distributions) == names(comb[j])) # nodo
            selected[[j]] = distributions[[node]][unlist(comb[j]),] # intervalo dentro del nodo
            names(selected)[j] = names(comb)[j]
          }
          # sel2 = unlist(selected, recursive = F)
          
          # combined factor's domain
          domain = getFactorDomain(selected)
          
          if(!is.null(domain)){
            # sel2$DomainFactorParents = domain
            # res[[k]] = sel2
            
            # Obtenemos la funciones de probabilidad
            arbol = lapply(selected, function(x){x[[length(x)]]})
            names(arbol) = names(selected)
            
            CPD_X_i = as.function(arbol[[paste0(X,"|",Y,":",className)]][[1]])
            CPD_Xi = as.function(arbol[[paste0(X,"|",className)]][[1]])
            integrando1 = function(x){
              # xx = seq(domain[1],domain[2],10^-4)
              # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
              return(CPD_X_i(x)*log(CPD_X_i(x)))
              # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
            }
            integrando2 = function(x){
              # xx = seq(domain[1],domain[2],10^-4)
              # minn = max(10^-5,min(CPD_X_i(xx),CPD_Xi(xx)))
              return(CPD_X_i(x)*log(CPD_Xi(x)))
              # return(CPD_X_i(x)*(log(CPD_X_i(x)/CPD_Xi(x))))
            }
            arbol$integrateX = tryCatch(integrate(integrando1,domain[1,X],domain[2,X])$value-
                                          integrate(integrando2,domain[1,X],domain[2,X])$value,
                                        error = function(e) {
                                          cat("Error:", e$message, "\n")
                                          0},warning = function(w) {
                                            cat("Advertencia:", w$message, "\n")
                                            0})
            arbol$integrateC = integrate.motbf(arbol[[paste0(className)]][[1]],
                                             domain[1,className],domain[2,className])
            
            MI = MI+arbol$integrateC*arbol$integrateX*arbol[[paste0(Y,"|",className)]]
            
          }
          
        }
      }
    }else{
      if(X_discrete){ # discreta|continua y clase continua
        selected = list()
        k = 1
        i = 1
        j = 1
        for(i in 1:nrow(combination.intervals)){
          comb = combination.intervals[i,, drop = FALSE]
          
          
          for(j in 1:length(comb)){
            node = which(names(distributions) == names(comb[j])) # nodo
            selected[[j]] = distributions[[node]][unlist(comb[j]),] # intervalo dentro del nodo
            names(selected)[j] = names(comb)[j]
          }
          # sel2 = unlist(selected, recursive = F)
          
          # combined factor's domain
          domain = getFactorDomain(selected)
          
          if(!is.null(domain)){
            # sel2$DomainFactorParents = domain
            # res[[k]] = sel2
            
            # Obtenemos la funciones de probabilidad
            arbol = lapply(selected, function(x){x[[length(x)]]})
            names(arbol) = names(selected)
            
            arbol$X = arbol[[paste0(X,"|",Y,":",className)]]*log(arbol[[paste0(X,"|",Y,":",className)]])-
              arbol[[paste0(X,"|",Y,":",className)]]*log(arbol[[paste0(X,"|",className)]])
            arbol$integrateC = integrate.motbf(arbol[[paste0(className)]][[1]],
                                             domain[1,className],domain[2,className])
            arbol$integrateY = integrate.motbf(arbol[[paste0(Y,"|",className)]][[1]],
                                             domain[1,Y],domain[2,Y])
            MI = MI+arbol$integrateC*arbol$X*arbol$integrateY
            
          }
          
        }
      }else{# continua|continua y clase continua
        # comprobar si cada combinacion de intervalos se puede multiplicar
        # y obtener el dominio sobre el cual la funcion es valida
        
        selected = list()
        k = 1
        i = 1
        j = 1
        for(i in 1:nrow(combination.intervals)){
          comb = combination.intervals[i,, drop = FALSE]
          
          
          for(j in 1:length(comb)){
            node = which(names(distributions) == names(comb[j])) # nodo
            selected[[j]] = distributions[[node]][unlist(comb[j]),] # intervalo dentro del nodo
            names(selected)[j] = names(comb)[j]
          }
          # combined factor's domain
          domain = getFactorDomain(selected)
          
          if(!is.null(domain)){
            arbol = lapply(selected, function(x){x[[length(x)]]})
            names(arbol) = names(selected)
            
            CPD_X_i = as.function(arbol[[paste0(X,"|",Y,":",className)]][[1]])
            CPD_Xi = as.function(arbol[[paste0(X,"|",className)]][[1]])
            integrando1 = function(x){
              return(CPD_X_i(x)*log(CPD_X_i(x)))
            }
            integrando2 = function(x){
              return(CPD_X_i(x)*log(CPD_Xi(x)))
            }
            arbol$integrateX = tryCatch(integrate(integrando1,domain[1,X],domain[2,X])$value-
                                          integrate(integrando2,domain[1,X],domain[2,X])$value,
                                        error = function(e) {
                                          cat("Error:", e$message, "\n")
                                          0},warning = function(w) {
                                            cat("Advertencia:", w$message, "\n")
                                            0})
            arbol$integrateC = integrate.motbf(arbol[[paste0(className)]][[1]],
                                             domain[1,className],domain[2,className])
            
            arbol$integrateY = integrate.motbf(arbol[[paste0(Y,"|",className)]][[1]],
                                             domain[1,Y],domain[2,Y])
            MI = MI+arbol$integrateC*arbol$integrateX*arbol$integrateY
            
          }
          
        }
      }
    }
  }
  return(list(MI=MI,distributions=distributions))
}

#' @noRd
# Funcion para calcular la informacion mutua de dos variables condicionada a otra (se exporta)--
cond_mut_information = function(data,className,varNames=NULL,
                                fit.args =NULL){
  
  # Comprobaciones
  # Chequeamos los argumentos para fijar las Mops------------
  fit.args = fit.args.null(fit.args)
  # browser()
  if(!(className%in%colnames(data))){
    stop(paste0("data must have ",className," as colnames",collapse = " "))
  }
  # C = data[,className]
  
  if(!is.null(varNames)){
    if(!all(varNames%in%colnames(data))){
      stop(paste0("data must have ",varNames," as colnames",collapse = " "))
    }
    data = data[,varNames]
  }
  
  if(length(data)!=3){
    stop("data must have 3 variables")
  }else{
    varNames = colnames(data)
  }
  varNames = setdiff(colnames(data),className)
  # Condiciona la segunda variable de varNames
  MI1 = cond_mut_information_compute(data[,c(varNames,className)],className,fit.args)
  # Condiciona la primera variable de varNames
  MI2 = cond_mut_information_compute(data[,c(varNames[2:1],className)],
                                     className,fit.args,
                                     distC = MI1$distributions[[className]],
                                     distY_C = MI1$distributions[[2]],
                                     distX_C = MI1$distributions[[4]])
  return((MI1$MI+MI2$MI)/2)
}



# Función para el algoritmo de chow-liu--------------------
## Función para el arbol maximal--------------------
#' @noRd
prim_maximal <- function(adj_mat,root = 1) {
  # Número de nodos en el grafo
  # browser()
  n <- nrow(adj_mat)
  if(is.null(n)){
    n=1# caso 1 variable
  }
  
  # Inicialización
  selected_nodes <- rep(FALSE, n)
  names(selected_nodes) = row.names(adj_mat)
  selected_nodes[root] <- TRUE  # Seleccionamos el nodo inicial arbitrario
  maximal_tree <- matrix(0, n, n,dimnames = dimnames(adj_mat))  # Para almacenar el árbol máximo
  
  # Bucle principal
  i = 1
  for (i in 1:(n-1)) {
    max_weight <- -Inf
    u <- -1#Fila
    v <- -1# Columna
    
    # Encontrar la arista de mayor peso que conecta un nodo seleccionado con uno no seleccionado
    for (j in which(selected_nodes)) {
      for (k in which(!selected_nodes)) {
        if (adj_mat[j, k] > max_weight) {
          max_weight <- adj_mat[j, k]
          u <- j
          v <- k
        }
      }
    }
    
    # Añadir la arista al árbol maximal
    if (u != -1 && v != -1) {
      maximal_tree[u, v] <- 1
      # maximal_tree[v, u] <- max_weight
      selected_nodes[v] <- TRUE
    }
  }
  return(maximal_tree)
}



# funcion para calcular la raiz de un TAN---------------------
#' @noRd
fitroot = function(MI){
  root = which.max(MI)%%nrow(MI)# Buscamos la fila que tiene el mayor valor de 
  # informacion mutua. Se obtiene como raiz la variable que condiciona en dicho valor
  root = ifelse(root==0,nrow(MI),root)# Al hacer el modulo, 0=root=nrow
  return(root)
}

## Example 1----------------------------------------------------

# set.seed(14)
# C = as.factor(sample(c("c1","c2"),100,T,c(0.4,0.6)))
# tbC = as.vector(table(C))
# set.seed(1452)
# X = c()
# X[C =="c1"] = rbeta(tbC[1],2,3)
# set.seed(548)
# X[C =="c2"] = 2+rbeta(tbC[2],2,3)
# set.seed(841)
# Y = as.factor(sample(c("y1","y2"),100,T))
# set.seed(74)
# Z = rnorm(100)
# data = data.frame(X,Y,Z,C)
# 
# target = "C"
# mutualInfoCond = NULL
# root=NULL
# distributions=NULL
# fit.args=NULL
# fit_tan(target = target,data = data,root=root)

## Example 2-----------------------------------------------------
# set.seed(14)
# C = as.factor(sample(c("c1","c2"),100,T,c(0.4,0.6)))
# tbC = as.vector(table(C))
# set.seed(1452)
# X = c()
# X[C =="c1"] = rbeta(tbC[1],2,3)
# set.seed(548)
# X[C =="c2"] = 2+rbeta(tbC[2],2,3)
# set.seed(841)
# Y = as.factor(sample(c("y1","y2"),100,T))
# set.seed(74)
# Z = rnorm(100)
# A= c()
# set.seed(12)
# A[C =="c1"] = sample(c("a1","a2"),tbC[1],T,c(0.7,0.3))
# set.seed(174)
# A[C =="c2"] = sample(c("a1","a2"),tbC[2],T,c(0.25,0.75))
# A = as.factor(A)
# target = "C"
# mutualInfoCond = NULL
# root=NULL
# distributions=NULL
# fit.args=fit.args.null(NULL)
# data = data.frame(X,Y,Z,C)
# data_new = data.frame(X,Y,Z,A,C)
# order = 1:4
# varsNew = "A"
# library(MoTBFs)
# library(logging)
# bn_old = fit_tan(target,data = data,all = T,parallel = F)
# mutualInfoCond_0 = bn_old$mutualInfoCond
# distributions_0 = bn_old$distributions
# 
# bnAll = fit_tan(target = target,data = data_new,root = root,
#                 distributions = distributions,
#                 mutualInfoCond = mutualInfoCond,
#                 all=T)
# dag=getDAG(bnAll$bn)
# graphviz.plot(dag)
# library(bnlearn)


#' @noRd
fit_nb = function(data,target,fit.args){
  # Seleccionamos las variables
  varNames = colnames(data)
  # variables predictoras
  varPred = varNames[varNames!=target]
  
  # Distribuciones
  distributions= fitNB_dis(data,target,fit.args)
  bn = getFormatedBN(distributions)
  return(bn)
}

#' @noRd
fitNB_dis = function(data,target,fit.args){
  # Seleccionamos las variables
  varNames = colnames(data)
  # variables predictoras
  varPred = varNames[varNames!=target]
  
  distributions = list()
  nDist = 1
  # Aprendizaje de la clase
  if(is.numeric(data[,target])){
    # PX <- univMoTBF(data[,Child], POTENTIAL_TYPE, maxParam=maxParam, scale = FALSE)
    PX = univMoTBF(x=data[,target],POTENTIAL_TYPE = fit.args$POTENTIAL_TYPE,
                   evalRange = fit.args$evalRange,nparam = fit.args$nparam,
                   maxParam = fit.args$maxParam,scale=fit.args$scale)
  }else{
    # Tan para clasificacion
    distC <- probDiscreteVariable(data[,target])
    distC <- list(Child=target, functions=list(distC), varType = "Discrete")
  }
  nDist = 1
  distributions[[nDist]] = distC
  names(distributions)[nDist]=paste0(target)
  ## Aprendizaje de las distribuciones X_i|C
  var_i = varPred[1]
  for(var_i in varPred){
    dist_i_C = do.call(conditionalMethod,
                       append(fit.args,
                              list(data=data,nameParents=target,
                                   nameChild=var_i)))
    nDist = nDist+1
    varType_i = ifelse(is.numeric(data[,var_i]),"Continuous","Discrete")
    dist_i_C = list(Child=var_i,functions=dist_i_C,varType = varType_i)
    distributions[[nDist]] = dist_i_C
    names(distributions)[nDist]=paste0(var_i,"|",target)
  }
  return(distributions)
}