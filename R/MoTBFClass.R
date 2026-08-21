### CLASS CONSTRUCTORS #####
# Low-level constructor for motbf objects
# 
# Defines an object of class \code{"motbf"}, os subclasses, and other basic functions for manipulating
# \code{"motbf"} objects.
# 
# @param x Preferably, a list containing an \code{'mte'} or \code{'mop'} univariate expression
# and other posibles elements like a \code{"numeric"} vector with the domain of the variable,
# the number of iterations needed to solve the optimization problem, among others.
# Any \R object can be entered, but the utility of this function is not to transform
# objects of other classes into objects of class \code{"motbf"}.
# @examples
# ## Subclass 'MOP'
# param <- c(1,2,3,4,5)
# MOPString <- asMOPString(param)
# fMOP <- motbf(MOPString)
# print(fMOP) ## fMOP
# as.character(fMOP)
# as.list(fMOP)
# is(fMOP)
# is.motbf(fMOP)
# 
# ## Subclass 'MTE'
# param <- c(6,7,8,9,10)
# MTEString <- asMTEString(param)
# fMTE <- motbf(MTEString)
# print(fMTE) ## MTE
# as.character(fMTE)
# as.list(fMTE)
# is(fMTE)
# is.motbf(fMTE)
#' @noRd
motbf <- function(x=0)
{
  if(!is.list(x)) x <- list(Function = noquote(x))
  # result <- x
  # class(result) <- c("univmotbf","motbf")
  # result
  structure(x, class = c("univmotbf", "motbf"))
}

#' @noRd
new_motbf <- function(x){
  if(!is.list(x)) x <- list(Function = noquote(x))
  structure(x, class = c("motbf"))
}

#' @noRd
new_univmotbf <- function(x){
  if(!is.list(x)) x <- list(Function = noquote(x))
  structure(x, class = c("univmotbf", "motbf"))
}

#' @noRd
new_mop <- function(x){
  if(!is.list(x)) x <- list(Function = noquote(x))
  structure(x, class = c("univmotbf", "mop", "motbf"))
}

#' @noRd
new_mte <- function(x){
  if(!is.list(x)) x <- list(Function = noquote(x))
  structure(x, class = c("univmotbf", "mte", "motbf"))
}

# Class \code{"piecewisemop"}
#' @noRd
new_piecewisemop <- function(x){
  structure(x, class = c("piecewisemop", "motbf"))
}

#' @noRd
new_jointmotbf <- function(x){
  if(!is.list(x)) x <- list(Function = noquote(x))
  structure(x, class = c("jointmotbf", "motbf"))
}


# Class \code{"motbf_fit"}
#' @noRd
new_motbf_fit <- function(x){
  structure(x, class = c("motbf_fit", "motbf"))
}

#' @noRd
new_motbf_fit_node <- function(x){
  structure(x, class = c("motbf.fit.node"))
}

# Class \code{"motbf_fit_cv"}
#' @noRd
new_motbf_fit_cv <- function(x){
  structure(x, class = c("motbf_fit_cv", "motbf"))
}



#' Class \code{"jointmotbf"}
#' 
#' DEPRECATED - Defines an object of class \code{"jointmotbf"} and other basic functions for 
#' manipulating \code{"jointmotbf"} objects.
#' 
# @name Class-JointMoTBF
# @rdname Class-JointMoTBF
#' @noRd
#' @param x Preferably, a list containing an expression
#' and other possible elements like a \code{"numeric"} matrix with the domain of the variables, 
#' the dimension of the variables, the number of iterations needed to solve the optimization problem,
#' among others. Any \R object can be entered, but the utility of this function is not to transform
#' objects of other classes into objects of class \code{"jointmotbf"}.
#' @seealso \link{jointMoTBF}
#' @examples
#' ## n.parameters is the product of the dimensions
#' dim <- c(3,3)
#' param <- seq(1,prod(dim), by=1)
#' ## Joint Function 
#' f <- list(Parameters=param, Dimensions=dim)
#' jointF <- jointMoTBF(f)
#' 
#' print(jointF) ## jointF
#' as.character(jointF)
#' is(jointF)
#' is.jointmotbf(jointF)
jointmotbf <- function(x = 0)
{
  if(!is.list(x)) x <- list(Function = noquote(x))
  result <- x
  class(result) <- c("jointmotbf", "motbf")
  result
}


## CLASS COERTION ########

#' Coerce MOTBF Objects to Character or Function
#'
#' Converts \code{'motbf'} and \code{'jointmotbf'} objects into character string expressions or executable R functions.
#'
#' @param x An object of class \code{'motbf'} or \code{'jointmotbf'}.
#' @param ... Further arguments passed to or from other methods. Not used currently.
#'
#' @return 
#' \itemize{
#'   \item \code{as.character}: Returns a character string representing the mathematical expression of the object.
#'   \item \code{as.function}: Returns an executable R function that accepts numeric arguments to evaluate the MoTBF expression.
#' }
#'
#' @name coercion-motbf
#' @rdname coercion-motbf
#' 
#' @examples
#' 
#' ## Example 1
#' X <- rchisq(5000, df = 3)
#' P <- univMoTBF(X, POTENTIAL_TYPE = "MOP"); P
#' as.function(P)(10)
#' 
#' ## Example 2
#' data <- data.frame(X = rnorm(100), Y = rexp(100))
#' dim <- c(3,2)
#' P <- jointmotbf.fit(data, dimensions = dim)
#' density <- as.function(P)(data[,1], data[,2])
#' sum(log(density))
#' 
#' @exportS3Method base::as.character motbf
as.character.motbf <- function(x, ...){ 
  as.character(x[[1]])
}

#' @rdname coercion-motbf
#' @exportS3Method base::as.function motbf
as.function.motbf <- function(x, ...){

  v <- getMotbfVar(x)

  f = eval(parse(text = paste("f <- function(",v,")",x)))
  formals(f) <- formals(f)[1:length(v)]
  names(formals(f)) <- v
  return(f)
  
}

#' @rdname coercion-motbf
#' @exportS3Method base::as.character jointmotbf
as.character.jointmotbf <- function(x, ...){
  as.character(x[[1]])
}


#' @rdname coercion-motbf
#' @exportS3Method base::as.function jointmotbf
as.function.jointmotbf <- function(x, ...)
{
  P <- x[[1]]
  # v <- getMotbfVar(P)
  v <- attr(x$Domain, "dimnames")[[2]]
  f= eval(parse(text = paste0("f <- function(",paste0(v, collapse = ', '),")",P)))
  
  formals(f) <- formals(f)[1:length(v)]
  names(formals(f)) <- v
  return(f)  
}




## CHECK CLASS ########

#' Check MoTBF Classes and Subclasses
#'
#' Utility functions to check whether an object belongs to a specific MoTBF class, or to identify its underlying subclass (\code{'mop'} or \code{'mte'}).
#'
#' @param x An object to be checked.
#' @param class Character string specifying the target class name.
#' @param fx An object to determine the subclass for.
#'
#' @return 
#' \itemize{
#'   \item \code{is.*}: Logical value (\code{TRUE} or \code{FALSE}) indicating if the object belongs to the checked class.
#'   \item \code{subclass}: A character string (\code{"mop"} or \code{"mte"}) specifying the underlying family of the object.
#' }
#'
#' @name is.motbf
#' @rdname is.motbf
#' @export
is.motbf <- function(x, class = "motbf"){ 
  is(x, class)
}

#' @rdname is.motbf
#' @export
is.univmotbf <- function(x, class = "univmotbf"){ 
  is(x, class)
}


#' @rdname is.motbf
#' @export
is.jointmotbf <- function(x, class = "jointmotbf"){ 
  is(x, class)
}

#' @rdname is.motbf
#' @export
is.motbf_fit <- function(x, class = "motbf_fit"){
  is(x, class)
}

#' @rdname is.motbf
#' @export
is.motbf_fit_cv<- function(x, class = "motbf_fit_cv"){
  is(x, class)
}

#' @rdname is.motbf
#' @export
is.mte <- function(x){
  if(is(x, 'mte')){
    return(TRUE)
  }
  subclass = tryCatch({x$Subclass}, error=function(e){NULL})
  if(!is.null(subclass)) return(subclass=="mte")
  f <- x[[1]] 
  l <- length(strsplit(as.character(f), split="exp", fixed=TRUE)[[1]])-1
  return(l!=0&&is.motbf(x))
}

#' @rdname is.motbf
#' @export
is.mop <- function(x){
  if(is(x, 'mop')){
    return(TRUE)
  }
  subclass = tryCatch({x$Subclass}, error=function(e){NULL})
  if(!is.null(subclass)) return(subclass=="mop")
  f <- x[[1]] 
  l <- length(strsplit(as.character(f), split="exp", fixed=TRUE)[[1]])-1
  return(l==0&&is.motbf(x))
}

#' @rdname is.motbf
#' @export
subclass <- function(fx)
{
  if(is.mop(fx)) return("mop")
  if(is.mte(fx)) return("mte")
}

