#' Joint MoTBF density learning
#' 
#' Function for learning joint MoTBFs.
#' The \code{jointmotbf.fit()} function is a wrapper of two internal (non-exported) functions: 
#' \code{getParamJoint()} and \code{fixParamJoint()}.
#' The first one gets the parameters by solving a quadratic optimization problem, minimizing
#' the mean squared error between the empirical joint CDF and the estimated CDF.
#' The density is obtained as the derivative of the estimated CDF.
#' The second one, \code{fixParamJoint()}, fixes the equation of the joint function using 
#' the previously learned parameters and converting this \code{"character"} string into an 
#' object of class \code{"jointmotbf"}.
#' 

#' @param X a dataset of class \code{"data.frame"}.
#' @param ranges a \code{"numeric"} matrix containing the range of the variables used to fit the function, 
#' where each column corresponds to a variable. If not specified, the range of each variable is computed from the data.
#' @param dimensions a \code{"numeric"} vector containing the number of parameters of each variable.
#' @param fitPoints an \code{"integer"} indicating the number of points per variable to use to build the expanded grid 
#' where the objective function will be evaluated when optimizing the parameters.
#' @param constraints an \code{"integer"} indicating the number of constraints under which to minimize the quadratic function.
#' @return
#' \code{jointmotbf.fit()} returns a list with the following elements: 
#' \item{\code{Function}}{The analytical expression of the learned density.}
#' \item{\code{Domain}}{A \code{"matrix"} containing the domain of each variable over which the density is defined.}
#' \item{\code{Iterations}}{The number of iterations needed to solve the problem.}
#' \item{\code{Time}}{The execution time.}
#' 
#' @importFrom Matrix nearPD
#' @examples
#
#' ## 1. EXAMPLE 
#' ## Generate a multinormal dataset
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100))
#' 
#' ## Joint learnings
#' dim <- c(2,3)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' 
#' P
#' attributes(P)
#' class(P)
#' 
#' ###############################################################################
#' ## MORE EXAMPLES ##############################################################
#' ###############################################################################
#' \donttest{
#' ## Generate a dataset
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100), X3 = rnorm(100))
#' 
#' ## Joint learnings
#' dim <- c(3,2,3)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' P
#' attributes(P)
#' class(P)
#' }
#' @export
#' 
jointmotbf.fit = function(X, ranges = NULL, dimensions = NULL, fitPoints = 10, constraints = 10){
  
  param <- getParamJoint(X, ranges, dimensions, fitPoints, constraints)
  P <- fixParamJoint(param)
  
  return(P)
}

#' @noRd
getParamJoint=function(X, ranges = NULL, dimensions = NULL, fitPoints = 10, constraints = 10){
  tm <- Sys.time()

  if(is.null(ranges)) ranges <- sapply(X, range)
  if(is.null(dimensions)){
    Fx <- lapply(X, univMoTBF, "MOP")
    ll <- list(); for(i in 1:length(Fx)) ll[[length(ll)+1]] <- length(coef(Fx[[i]]))
  } else {
    ll <- list(); for(i in 1:length(dimensions)) ll[[length(ll)+1]] <- dimensions[i]
  }
  
  dim <- c(); pos <- which(c(length(ll[[1]]), length(ll[[2]]))==min(length(ll[[1]]), length(ll[[2]])))
  if(pos[1]==1) {pos1 <- 1; pos2 <- 2}
  if(pos[1]!=1) {pos1 <- 2; pos2 <- 1}
  for(i in 1:length(ll[[pos1]])) for(j in 1:length(ll[[pos2]])) dim <- c(dim,list(c(ll[[pos1]][i],ll[[pos2]][j])))
  if(pos1==2) for(i in 1:length(dim)) dim[[i]] <- dim[[i]][length(dim[[i]]):1]
  
  size <- prod(sapply(1:length(ll), function(i) length(ll[[i]])))
  si <- lapply(1:ncol(X),function(i) rep(ll[[i]][1], size))
  dim <- c(); for(i in 1:length(si)) dim <- cbind(dim,unlist(si[[i]]))
  dim <- lapply(1:nrow(dim), function(i) dim[i,])
  
  ## Create the grid: depending on the n.record and the n.variables
  eg <- list()
  for(i in 1: ncol(X)){
    eg[[i]] <- seq(ranges[1,i], ranges[2,i], length.out = fitPoints)
    
  }
  x <- expand.grid(eg)
  
  ## Cumulative densities
  y <- jointCDF(X,x)
  
  P <- c()
  for(s in 1:length(dim)){
    Xt <- ranges; pos <- coefExpJointCDF(dim[[s]])
    n <- nrow(x); nterms <- length(pos)
    
    
    x1 <- list()
    for(j in 1:ncol(X)){
      pos2 <- c(0,unlist(lapply(pos, function(x){x[j]})))
      x1[[j]] <- outer(x[,j], pos2, FUN = "^")
    }
    xx <-Reduce("*", x1)
    
    ## Dmat
    XX <- t(xx)%*%xx 
    
    ## dvec
    Xy <- t(xx)%*%y 
    
    ## Constrains
    xNew <- list()
    for( j in 1:ncol(Xt)){
      xnew <- seq(min(Xt[,j]), max(Xt[,j]), length=constraints)
      if(xnew[length(xnew)]!=max(Xt[,j])) xnew <- c(xnew,max(Xt[,j]))
      xNew[[length(xNew)+1]] <- xnew
    }
    xN <- expand.grid(xNew)
    
    ma22 <- c()
    for(i in 1:ncol(X)){
      d <- sapply(1:nterms, function(j) pos[[j]][i]*(xN[,i]^(pos[[j]][i]-1)))
      tf <- sapply(1:length(pos), function(j) pos[[j]][i]==0)
      if(any(tf)) d[,which(tf)] <- 0
      d <- cbind(rep(0, nrow(xN)),d)
      ma22[[length(ma22)+1]] <- t(d)
    }
    ma2 <- Reduce("*", ma22)
    
    ma3 <- c()
    for(j in 1:ncol(Xt)){
      ma <- 0
      for(i in 1:nterms) ma <- rbind(ma,(max(Xt[,j])^pos[[i]][j]) - (min(Xt[,j])^pos[[i]][j]))
      ma3 <- cbind(ma3,ma)
    }
    for(i in 1:(ncol(X)-1)) ma3[,i+1] <- ma3[,i] * ma3[,i+1]
    ma33 <- cbind(ma3[,i+1])
    
    ## Amat
    AA <- cbind(ma33, ma2)
    
    ## bvec
    B <- c(1,rep(1.0E-5,ncol(ma2)))
    
    tr <- tryCatch(solve.QP(XX, Xy, AA, B, meq=1), error = function(e) NULL)
    if(is.null(tr)==TRUE){
      message("matrix D in quadratic function has been approximated to the nearest positive definite!")
      tr <- solve.QP(nearPD(XX)$mat, Xy, AA, B, meq=1)
      # return(NULL)
    }   
    
    finaltm <- Sys.time() - tm
    soluc <- tr 
    parameters <- soluc$solution
    
    # if(is.null(tr)==TRUE){
    #   return(NULL)
    # }else{    
    #   soluc <- tr 
    #   parameters <- soluc$solution
    # }
    
    ## Parameters PDF
    a <- sapply(pos, prod)
    param <- parameters[-1]*a
    parameters <- param[which(param!=0)]
    
    P <- list(Parameters = parameters, 
              Dimension = dim[[s]], Range = Xt, 
              Iterations = tr$iterations[1], 
              Time=finaltm)
  }
  
  return(P)
}


#' @noRd
fixParamJoint <- function(object){
  parameters <- object$Parameters
  dimensions <- object$Dimension

  v = colnames(object$Range)
  if(length(parameters)==1) {
    s <- paste(parameters[1], "+0", sep="")
    for(i in 1:length(v)) s <- paste(s, "*", v[i], sep="")
    P <- list(Function = s, Domain = object$Range)
    P <- new_jointmotbf(P)
    return(P)
  }
  pos <- coefExpJointCDF(dimensions)
  t <- list()
  for(i in 1:length(pos)) if(prod(pos[[i]])!=0) 
    t[[length(t)+1]] <- pos[[i]] else next
  pos <- t
  
  if(any(dimensions==1)) posvardim1 <- which(dimensions==1)
  
  if(parameters[1]==0) s <- c() else s <- parameters[1]
  for(i in 2:length(parameters)){
    if(is.null(s)){
      s <- paste(s, parameters[i],sep="")
    } else{
      s <- paste(s, ifelse(parameters[i]>=0, "+", ""), parameters[i],sep="")
      if(i != length(parameters)) {
        for(p in 1:length(v)){
          s <- paste(s, ifelse((pos[[i]][p]-1)==0, "", paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "", paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
        } 
        
      }
      if(i == length(parameters)){
        for(p in 1:length(v)) s <- paste(s, ifelse((pos[[i]][p]-1)==0, ifelse(any(p==posvardim1),
                                                                              paste("*",v[p],"^",pos[[i]][p]-1, sep=""),""), 
                                                   paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "",
                                                                          paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
      }
    }
  }
  P <- noquote(s)
  rango = object$Range
  attr(rango, "modelVars") = v
  P <- list(Function = P, Domain = rango,
            Iterations = object$Iterations, Time = object$Time)
  P <- new_jointmotbf(P)
  return(P)
}


#' Degree Function
#'
#' Compute the degree for each term of a joint CDF.
#'
#' @param dimensions A \code{"numeric"} vector including the number of parameters of each variable.
#' @return A list with n element. Each element contains a \code{numeric} vector with the degree for 
#' each variable and each term of the joint CDF.
#' @noRd
#' @examples
#'
#' ## Dimension of the joint PDF of 2 variables
#' dim <- c(4,5) 
#' ## Potentials of each term of the CDF
#' d <- coefExpJointCDF(dim)
#' length(d) + 1 ## plus 1 because of the constant coefficient
#'
#' ## Dimension of the joint density function of 2 variables
#' dim <- c(5,5,3)
#' ## Potentials of the cumulative function
#' coefExpJointCDF(dim)
#'
coefExpJointCDF <- function(dimensions){
  t <- lapply(1:length(dimensions), function(i) 0:dimensions[i])
  l <- sapply(1:length(t), function(i) length(t[[i]]))
  mm <- c()
  for(i in 1:length(t)){
    m <- c()
    for(j in 1:length(t[[i]])) m <- c(m,rep(t[[i]][j], prod(l[-c(1:i)])))
    if(i!=1) m <- rep(m, l[i-1])
    mm <- cbind(mm, m)
  }
  mm <- mm[-1,]
  colnames(mm) <- NULL
  if(is.matrix(mm)) pos <- lapply(1:nrow(mm), function(i) mm[i,])
  else pos <-list(as.numeric(mm))
  return(pos)
}


#' Joint MoTBFs CDFs
#' 
#' Function to compute multivariate Cumulative Distribution Functions. 
#' 
#' @param df The dataset as an object of class \code{data.frame}.
#' @param grid a \code{data.frame} with the selected data points where the objective function
#' will be evaluated when optimizing the parameters.
#' @return \code{jointCDF()} returns a vector.
#' @examples
#' 
#' ## Create dataset with 2 variables
#' n = 2
#' size = 50
#' df <- as.data.frame(matrix(round(rnorm(size*n),2), ncol = n))
#' 
#' ## Create grid dataset
#' npointsgrid <- 10
#' ranges <- sapply(df, range)
#' eg <- list()
#' for(i in 1: ncol(df)){
#'   eg[[i]] <- seq(ranges[1,i], ranges[2,i], length.out = npointsgrid)
#' }
#' 
#' x <- expand.grid(eg)
#' 
#' ## Joint cumulative values
#' jointCDF(df = df, grid = x)
#' 
#' @noRd
jointCDF <- function(df, grid){
  # apply(grid, MARGIN=1, posGrid_fast, df=df)
  x <- grid
  n <- ncol(df)
  if(ncol(df)>10){
    apply(x, MARGIN=1, 
          FUN = function(x,df){
            b = which(colSums(t(df)<=x)==n)
            out = length(b)/nrow(df)
          }, 
          df=df)
  }else{
    apply(grid, MARGIN=1, 
          FUN = function(grid,df){
            n <- ncol(df)
            b = switch(n-1,
                       which(df[,1]<=grid[1] & df[,2]<=grid[2]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2]&df[,3]<=grid[3]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8] & df[,9]<=grid[9]),
                       which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8] & df[,9]<=grid[9] & df[,10]<=grid[10])
            )
            out = length(b)/nrow(df)
          }, 
          df=df)
  }
} 


# Alternative to jointCDF, slower.
# jointCDF <- function(df, grid){
#   n <- ncol(df)
#   if(ncol(df)>10){
#     apply(grid, MARGIN=1, 
#           FUN = function(grid,df){
#             b = which(colSums(t(df)<=grid)==n)
#             out = length(b)/nrow(df)
#           }, 
#           df=df)
#   }else{
#     apply(grid, MARGIN=1, 
#           FUN = function(grid,df){
#             n <- ncol(df)
#             b = switch(n-1,
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2]&df[,3]<=grid[3]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8] & df[,9]<=grid[9]),
#                        which(df[,1]<=grid[1] & df[,2]<=grid[2] & df[,3]<=grid[3] & df[,4]<=grid[4] & df[,5]<=grid[5] & df[,6]<=grid[6] & df[,7]<=grid[7] & df[,8]<=grid[8] & df[,9]<=grid[9] & df[,10]<=grid[10])
#             )
#             out = length(b)/nrow(df)
#           }, 
#           df=df)
#   }
# } 
# jointCDF <- function(data, grid) unlist(lapply(1:nrow(grid), posGrid, data, grid))
# 
# posGrid <- function(nrows, data, grid)
# {
#   p <- lapply(1:ncol(grid), function(i) which(data[,i]<=grid[nrows,i]))
#   pp <- c()
#   for(i in 1:length(p)){
#     pp <- c(pp,p[[i]])
#     if(i!=1) pp <- pp[duplicated(pp, last=TRUE)]
#   }
#   return(length(pp)/nrow(data))
# }



#' Coefficients of a \code{"jointmotbf"} object
#' 
#' Extracts the parameters of a joint MoTBF density. 
#' 
#' @param object An object of class \code{"jointmotbf"}.
#' @param \dots Other arguments, unnecessary for this function.
#' @return A \code{"numeric"} vector with the parameters of the function.
#' @seealso \link{jointmotbf.fit} 
#' @export
#' @examples
#' ## Generate a dataset
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100))
#' 
#' ## Joint function
#' dim <-c(2,4)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' P$Time
#' 
#' ## Coefficients
#' coef(P)
#'

coef.jointmotbf <- function(object, ...){ 
  coeffMOP(object)
}




#' Extract Variables of MoTBFs
#' 
#' Get the names of the variables involved in a  \code{"univmotbf"} or \code{jointmotbf} object.
#' 
#' @param P An object of class \code{"univmotbf"} or \code{"jointmotbf"}.
#' @return A \code{"character"} vector with the names of the variables in the function.
#' @export
#' @examples
#' 
#' # 1. EXAMPLE
#' ## Generate a dataset
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100))
#' 
#' ## Joint function
#' dim <-c(3,2)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' P
#' 
#' ## Variables
#' getMotbfVar(P)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' \donttest{
#' ## Generate a dataset
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100), X3 = rnorm(100))
#' 
#' ## Joint function
#' dim <- c(2,1,3)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' 
#' ## Variables
#' getMotbfVar(P)
#' }

getMotbfVar <- function(P){
  

  if(is.jointmotbf(P)||is.mop(P)){
    stmop = splitMOP(P) # mop split into coefficients and literal part
  }else if(is.mte(P)){
    stmop = splitMTE(P)
  }else{
    stop("Fail to detect motbf class")
  }
  
  res = stmop$Variables
  
  return(res)
}

#' Extract Dimension of MoTBFs
#' 
#' Get the dimension of \code{"motbf"} and \code{"jointmotbf"} densities.
#' 
#' @param P An object of class \code{"motbf"} and subclass 'mop' or \code{"jointmotbf"}.
#' @return Dimension of the function.
#' @seealso \link{univMoTBF} and \link{jointmotbf.fit}
#' @export
#' @examples
#' ## 1. EXAMPLE 
#' ## Data
#' X <- rnorm(2000)
#' 
#' ## Univariate function
#' f <- univMoTBF(X, POTENTIAL_TYPE = "MOP")
#' getMotbfDim(f)
#' 
#' ## 2. EXAMPLE 
#' ## Dataset with 2 variables
#' X <- data.frame(x = rnorm(100), y = rnorm(100))
#' 
#' ## Joint function
#' dim <- c(2,3)
#' P <- jointmotbf.fit(X, dimensions = dim)
#' 
#' ## Dimension of the joint function
#' getMotbfDim(P)
#'

getMotbfDim <- function(P){
  
  if(is.jointmotbf(P)||is.mop(P)){
    stmop = splitMOP(P) # mop split into coefficients and literal part
  }else if(is.mte(P)){
    stmop = splitMTE(P)
  }else{
    stop("Fail to detect motbf class")
  }
  
  
  b = stmop$Exponents
  
  b2 = strsplit(b, split = '\\*')
  
  b3 = unique(unlist(b2))
  
  b4 = strsplit(b3, split = '\\^')
  
  b5 = data.frame(var = sapply(b4, '[',1), grade = ifelse(is.na(sapply(b4, '[',2)),1,sapply(b4, '[',2)))
  
  b6 = b5[which(!is.na(b5$var)),]
  
  b7 = stats::aggregate(b6, by = list(b6$var), FUN = max)
  
  dim = as.numeric(b7$grade)+1
  names(dim) = b7$var
  
  return(dim)
  
}



#' Evaluation of MoTBFs
#' 
#' Evaluates a univariate (\code{"mop"}, \code{"mte"}) or joint distribution (\code{"jointmotbf"}) at a specific point.
#' 
#' @param P An object of class \code{"univmotbf"} or \code{"jointmotbf"}.
#' @param values A list with the name of the variables equal to the values to be evaluated.
#' @return If all the variables in the equation are evaluated, then a \code{"numeric"} value
#' is returned. Otherwise, an \code{"univmotbf"} object or a \code{"jointmotbf"} object is returned.
#' @export
#' @examples
#' #' ## 1. EXAMPLE
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rexp(100))
#' 
#' ## Joint function
#' dim <- c(3,3) # dim <- c(5,4)
#' P <- jointmotbf.fit(X, dimensions = dim)
#' 
#' ## Evaluation
#' val <- list(x = -1.5, y = 3)
#' eval.motbf(P, values = val)
#' val <- list(x = -1.5)
#' eval.motbf(P, values = val)
#' val <- list(y = 3)
#' eval.motbf(P, values = val)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ############################################################################## 
#' \donttest{
#' ## Dataset with 3 variables
#' X <- data.frame(x = rnorm(100), y = rexp(100), z = rnorm(100, 1))
#' 
#' ## Joint function
#' dim <- c(2,1,3)
#' P <- jointmotbf.fit(X, dimensions = dim)
#' P
#' 
#' ## Evaluation
#' val <- list(x = 0.8, y = -2.1, z = 1.2)
#' eval.motbf(P, values = val)
#' val <- list(x = 0.8, y = 1.2)
#' eval.motbf(P, values = val)
#' val <- list(y = -2.1)
#' eval.motbf(P, values = val)
#' }


eval.motbf = function(P, values){
  # Check input arguments
  if(!is.motbf(P)){
    stop('Argument "P" is not of class "motbf"')
  }
  if(is.univmotbf(P)){
    f = as.function(P)
    values = unname(unlist(values))
    return(f(values))
  }
  if(!is.list(values)||is.null(names(values))){
    stop('Argument "values" is not a named list')
  }
  
  # Get coefficients and literal part of polynomial
  A = splitMOP(P)
  coeffmop = A$Coefficients
  expo = A$Exponents
  
  lp = strsplit(expo, split="(\\*)", perl = TRUE) #literal part
  var <- names(values)
  
  for(j in 1:length(var)){# loop over observed variables in argument 'values'
    x = var[j] # character name of variable
    
    # extract exponent of variable 'x' in each term (and store it in 'b')
    tm1 = sapply(lp, grep, pattern = x, value = TRUE) 
    
    b = unname(sapply(tm1, function(x){if(length(x)==0){""}else{x}}))# if variable is not in term, leave it blank
    b = ifelse(b == "", 0, b);b # exponent for independent term (substitute "" for 0)
    b = gsub(paste0(x, '|\\^'),'',b) # get exponents
    b = as.numeric(ifelse(b == '', 1, b)) # when exponent is 1, it is shown as "" (substitute "" for 1)
    
    val = values[[x]] # value to substitute variable for
    coeffmop = coeffmop*val^b
  }
  
  mopvariables = A$Variables # all variables in the joint distribution
  
  if(all(mopvariables%in% var)){ # if function is evaluated in all variables, the result is numeric 
    return(sum(coeffmop))
    
  }else{# Non-observed variables
    nonObservedVars  = mopvariables[which(!(mopvariables %in% var))]
    litP = c()
    for(i in 1:length(nonObservedVars)){# loop over Non-observed variables to get the literal part of the evaluated polynomial
      x = nonObservedVars[i]
      tm1 = sapply(lp, grep, pattern = x, value = TRUE);tm1
      b = unname(sapply(tm1, function(x){if(length(x)==0){""}else{x}}));b
      
      LPi = ifelse(b =='', '', paste0("*",b));LPi # literal part in iteration i
      
      litP = paste0(litP, LPi)
    }
    # Sum terms with equal literal part
    terms2sum = unique(litP)
    
    coefSum = c()
    for(s in 1:length(terms2sum)){
      coefSum[s] = sum(coeffmop[which(litP == terms2sum[s])])
    }
    sn = ifelse(coefSum<0, "", '+')
    sn[1] = ''
    
    res = paste0(sn, coefSum, terms2sum, collapse = "")
    
    # Return result
    if(length(nonObservedVars)>1){
      f = new_jointmotbf(res)
    }else{
      f = new_mop(res)
    }
    return(f)
  }
}



#' Marginalization of MoTBFs
#'
#' Computes the marginal densities from a \code{"jointmotbf"} object.
#'
#'@param P An object of class \code{"jointmotbf"}, i.e., the joint density function.
#'@param var A vector containing the names of the marginal variables (those being retained in 'P').
#' This argument accepts the \code{"numeric"} position (w.r.t. P$Domain) or the \code{"character"} name of the variables.
#'@return 
#'If 'var' is an atomic vector, the marginal distribution of the variable specified in 'var' is computed 
#'from the joint distribution and the result is an object of class \code{"univmotbf"}, i.e., a univariate distribution. 
#'
#'If 'var' contains more than one variable, the joint distribution over this subset is computed, i.e., 
#'the variables NOT contained in 'var' are marginalized out and the result is an object of class \code{"jointmotbf"}.
#'@seealso \link{jointmotbf.fit} and \link{eval.motbf}
#'@export
#'@examples
#' ## 1. EXAMPLE 
#' ## Dataset with 2 variables
#' data <- data.frame(x = rnorm(100), y = rnorm(100))
#' 
#' ## Joint function
#' dim <- c(4,3)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' 
#' ## Marginal
#' marginal.jointmotbf(P, var = "x")
#' marginal.jointmotbf(P, var = 2)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' \donttest{
#' ## Generate a dataset with 3 variables
#' data <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#' 
#' ## Joint function
#' dim <- c(2,2,3)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' 
#' ## Marginal
#' marginal.jointmotbf(P, var = "x")
#' marginal.jointmotbf(P, var = "y")
#' marginal.jointmotbf(P, var = c("x", "z"))
#' }

marginal.jointmotbf <- function(P, var){
  
  suppressWarnings({
    varN <- colnames(P$Domain) # todas las variables del MOP
    if(is.numeric(var)){
      var = colnames(P$Domain)[var]
    }
    noVar=varN[!varN%in%var] # variables que se marginalizan
    l <- length(noVar)
    if(l == 0){
      return(P)
    }
    
    for(j in 1:l){
      varN <- colnames(P$Domain) # todas las variables del modelo
      posnoVar <- which(varN%in%noVar[j])
      posvar <- which(varN%in%var)
      
      min <- lapply(posnoVar, function(i) P$Domain[1,i])
      names(min) <- noVar[j]
      max <- lapply(posnoVar, function(i) P$Domain[2,i])
      names(max) <- noVar[j]
      
      iP <- integrateJointmotbf(P, noVar[j])
      minF <- eval.motbf(iP, min)
      maxF <- eval.motbf(iP, max)
      if(is.numeric(minF)&&is.numeric(maxF)){
        P <- list(Function=(maxF-minF), Subclass="mop", Domain= P$Domain[,which(!varN%in%noVar[j]), drop = FALSE])
        P <- new_mop(P)
        return(P)
      }
      
      #----#
      ## Min part
      minPart = splitMOP(minF)
      
      
      ## Max part
      maxPart = splitMOP(maxF)
      
      
      if(!all(minPart$Exponents == maxPart$Exponents)){
        stop('Some problem ocurred marginalizing')
      }
      newcoef = maxPart$Coefficients - minPart$Coefficients
      
      st = paste0(newcoef[1],maxPart$Exponents[1],paste0(ifelse(sign(newcoef[-1])<0,"", "+"), newcoef[-1], maxPart$Exponents[-1], collapse = ""))
      

      if(length(getMotbfVar(new_motbf(st)))>1){
        rango = P$Domain[,-which(colnames(P$Domain)%in%noVar[j]), drop = FALSE]
        P <- list(Function=noquote(st), Domain= rango)
        P <- new_jointmotbf(P)
      }else{
        rango = P$Domain[,-which(colnames(P$Domain)%in%noVar), drop = FALSE]
        P <- list(Function=noquote(st), Subclass="mop", Domain= rango)
        P <- new_mop(P)
        if(length(unique(coeffPol(P)))==1){
          st <- paste(sum(coef(P)), ifelse(unique(coeffPol(P))==0, "",paste("*", var, ifelse(unique(coeffPol(P))==1, "", paste("^", unique(coeffPol(P)), sep="")), sep="")), sep="")
          P <- list(Function=noquote(st), Subclass="mop", Domain= P$Domain)
          P <- new_mop(P)
        }
      }
    }# end loop over j (variables to marginalize out)
    
    return(P)
  })
}

