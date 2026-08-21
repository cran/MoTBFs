#' Variable selection for MoTBFs
#' 
#' Perform variable selection for Bayesian networks of class MoTBF.
#' 
#' @param data  an object of class \code{"data.frame"}, which can contain continuous and discrete variables.
#' @param dag a character string indicating the structural learning algorithm to be applied to the training data. 
#' Available options are naive Bayes (NB), tree augmented naive Bayes (TAN) and hill-climbing (HC).
#' @param loss a character string indicating which loss function should be used. 
#' Currently, two options are available: 'logl', for the log-likelihood of the model; and 'pred', for the predictive error. See details.
#' @param target a character string indicating which node is the target. 
#' @param order a character vector indicating the order in which the predictive variables are included in the model. If
#' it is NULL, the predictors are ordered regarding their mutual information with the target variable.
#' @param method a character vector indicating the method used for the variable selection. Available methods are 'forward' (default),
#' 'gf' (greedy forward selection) and 'iwsr' (Incremental Wrapper Sequential Subset with Replacement).
#' @param fit.args a list containing optional arguments used to fit the models. These arguments must be those accepted by function 
#' \link{motbf.fit}, i.e., 'numIntervals' (4), 'POTENTIAL_TYPE' ('MOP'), 'maxParam' (NULL), 's' (NULL), 'priorData' (NULL) or 'scale' (TRUE). 
#' If fit.args is left NULL, the default values (in brackets) for those arguments will be used.
#' @param loss.args a list containing optional arguments related to the loss functions. Currently available arguments are:
#'  'loss.matrix' and 'percentage_test'. See details.
#' @param cv.args a list containing optional arguments related to the cross-validation method. Currently available arguments are:
#' 'k', 'seed'. See details.
#' @param verbose Logical; if \code{TRUE}, prints execution messages and progress updates to the console.
#' @details 
#' A filter-wrapper variable selection procedure is implemented. 
#' Firstly, the explanatory variables are ordered according to their mutual information
#' with the target variable, unless argument \code{'order'} is not null, in which case
#' the given order is followed. 
#' Then, the first variable in the ordered set and the target are used to fit an initial model.
#' The model is validated by means of a k-fold cross validation (`cv.args[[k]] >= 2`). 
#' However, it is possible to validate on a test set (`cv.args[[k]] = 1`), or on the same training set (`cv.args[[k]] = 0`).
#' The model is evaluated in terms of its log-likelihood (loss = 'logl') or its predictive accuracy (loss = 'pred').
#' In the latter case, the root mean squared error is computed for regression models, 
#' whereas the classification accuracy is computed for classification models.
#' Afterwards, the remaining explanatory variables are included in the model, one by one, according to the aforementioned order,
#' and a new model is obtained. Whenever the inclusion of a variable increases the accuracy of the model 
#' (increases its log-likelihood or classification error, or decreases its root mean squared error), it is kept; otherwise,
#' it is excluded from the model.
#' 
#' 
#' Details on the \code{loss} argument:
#'\describe{ 
#' \item{\code{'logl'}}{The log-likelihood of the model is computed. This measure is available for both basis functions, MTE and MOP.}
#' \item{\code{'pred'}}{This option is only available for MOPs. 
#' The predictive error (root mean squared error) or classification accuracy is computed as the loss function, 
#' depending on the nature of the target variable (continuous or discrete, respectively).
#' The program will guess which type the target variable is and compute the corresponding measure.}
#' }
#' 
#' Details on the \code{loss.args} argument. Currently, two arguments can be specified within this list:
#' \describe{
#' \item{\code{'p.test'}}{Only used if \code{k = 1}. 
#' This argument specifies the proportion of the data set that goes to the test set (between 0 and 1).}
#' \item{\code{'loss.matrix'}}{A squared \code{matrix} used to compute a weighted classification accuracy for discrete targets. 
#' This matrix is multiplied by the confusion matrix element by element, 
#' and the resulting matrix is used to compute the classification accuracy. 
#' The loss matrix allows to increase the penalty of some user-specified errors. 
#' If the loss matrix provided is a constant matrix of ones, the result is the standard classification accuracy.
#' Note that the diagonal of the loss matrix is regarded as a reward, while the off-diagonal is regarded as a cost.}
#' }
#' 
#' Details in the \code{cv.args} argument. Currently, two arguments can be specified within this list:
#' \describe{
#' \item{\code{k}}{an integer indicating the number of folds to split the data set. If k =  0, the train and test sets
#' are the same data; if k = 1, hold-out validation is carried out, i.e., the data set is split in train (80% by default) and test (20% by default);
#' finally, if k >=2, k-fold cross validation is carried out.}
#' \item{\code{seed}}{an integer to specify the seed. The k-folds are created randomly, so one might expect slightly different results unless 'seed' is used.}
#' }
#' 
#' @return A list of 4 elements
#' \describe{
#' \item{data}{a \code{"data.frame"} of the selected variables, including the target.}
#' \item{index}{the index of the selected variables (with respect to the input dataset), including the target.}
#' \item{loss}{the loss value of the best model.}
#' \item{crossvalidation}{an object of class \code{"motbf.fit.cv"}, containing the result of the cross-validation of the best model.}
#' }
#' @export
#' @examples 
#' #################
#' ### EXAMPLE 1 ###
#' #################
#' # Perform variable selection on the iris dataset.
#' # Use the TAN structure as DAG and the Species variable as target.
#' vs = variableSelection(iris, dag = 'TAN', loss = 'pred', target = 'Species', 
#'   cv.args = list(k = 1, seed = 1023))
#' 
#' # Check out the results of the best model.
#' summary(vs$crossvalidation)
#' 
#' 
variableSelection <- function(data, dag, loss, target = NULL, order = NULL, 
                              method = 'forward',
                              fit.args = NULL, loss.args = NULL, cv.args = NULL, verbose = TRUE){
  # check input arguments
  if(!(dag %in% c('NB', 'TAN', 'HC'))){
    stop("Argument 'dag' must be one of the following options: 
         'NB' (naive Bayes structure); 
         'TAN' (tree augmented Naive Bayes structure); 
         'HC' (structure via hill-climbing algorithm).")
  }
  
  if(is.null(target)){
    stop("Argument 'target' is empty with no default.")
  }
  # if(!((is.null(target) & dag == 'HC') & (is.null(target) & loss == 'logl'))){
  #   stop("Argument 'target' is empty with no default. Argument 'target' is mandatory when dag = 'NB' or 'TAN', or when loss = 'pred'. ")
  # }
  
  if(!(target%in%colnames(data))){
    stop("Argument 'target' is not available in the dataset 'data'. ")
  }
  
  if(!(loss %in% c('logl', 'pred'))){
    stop("Argument 'loss' must be one of the following options: 
         'logl' (loglikelihood of the model); 
         'pred' (predictive error).")
  }
  
  if(!(method %in% c("forward", "gf","iwsr"))){
    stop(paste0("Selective method '", method, "' not found. Available options are: 'forward', 'gf' (greedy forward) and 'iwsr' (Incremental Wrapper Sequential Subset with Replacement ."))
    
  }
  
  target_index = which(colnames(data) == target)
  
  # ORDER EXPLANATORY VARIABLES, according to their mutual information with the class (if order = NULL)
  if(is.null(order)){
    # Discretize the variables
    datD = discretizeVariablesEWdis(data, 5, factor = TRUE)
    # datD = infotheo::discretize(data)
    mi = infotheo::mutinformation(datD);mi
    mi = mi[,target]
    
    order = sort(mi[-target_index], decreasing = TRUE)
    order = match(names(order),colnames(data))
  }else{
    if(!all(order %in% colnames(data))){
      stop("Some explatatory variables included in the 'order' vector are not present in the data.frame 'data'. ")
    }
    order = match(order,colnames(data))
  }
  
  if(!is.null(cv.args$k)){
    k = cv.args$k
  }else{
    k = 10
  }
  if(!is.null(cv.args$seed)){
    seed = cv.args$seed
  }else{
    seed = NULL
  }
  
  # STORE CROSS-VALIDATION ARGUMENTS
  args = list(dag = dag, loss = loss, k = k, target = target, seed = seed, fit.args = fit.args, loss.args = loss.args)
  
  # BUILD INICIAL MODEL
  
  # Subset data.frame
  varsModelo = order[1]
  dfModelo <- data[,c(varsModelo, target_index)]
  
  # Perform k-fold cross validation
  motbf.cv.args = append(list(data = dfModelo), args)
  cv.mse = do.call(motbf.cv, motbf.cv.args)
  # cv.mse = motbf.cv(data = dfModelo, dag = dag, loss = loss, k = k, target = target, seed = seed, fit.args = fit.args, loss.args = loss.args)
  
  # Model accuracy
  as_actual = mean(sapply(cv.mse, '[[', 'loss'))
  
  # DETECT WHETHER TO MAXIMIZE OR MINIMIZE
  if(loss == 'logl'){
    loss_info = 'Log-likelihood'
    # the higher, the better
    maximize = TRUE
  }else if(loss == 'pred'){
    if(attr(cv.mse, 'target.type') == 'Discrete'){
      loss_info = 'Classification accuracy'
      maximize = TRUE
    }else{
      loss_info = 'Root mean squared error'
      # the lower, the better
      maximize = FALSE
    }
  }
  
  if(verbose){
    cat("Initial model:", paste0(colnames(dfModelo), collapse = ', '), "\n ", loss_info ,"=", as_actual, "\n\n")
  }
  
  globalImprove = TRUE
  while(globalImprove){
    
    globalImprove = FALSE
    for(i in 2:length(order)) {# Introducir 1 a 1 cada variable y calcular el accuracy.
    
    bestSet = NULL
    
    newVar = order[i] 
    
    if(newVar %in% varsModelo){
      next # saltar variable que ya está en el modelo
    }
    
    if(verbose){
      cat('NEW variable:', colnames(data)[newVar], '\n')
    }
    
    # REPLACEMENT
    ##############-
    if(method == 'iwsr'){
      if(verbose){
        cat("REPLACEMENT \n")
      }
      for(j in 1:length(varsModelo)){
        varsModeloNew = varsModelo
        if(verbose){
          cat("Replace", colnames(data)[varsModeloNew[j]], '\n')
        }
        varsModeloNew[j] = newVar
        
        dfModelo = data[,c(varsModeloNew, target_index)]
        
        # FIT NEW MODEL
        motbf.cv.args = append(list(data = dfModelo), args)
        cv.mse_new = do.call(motbf.cv, motbf.cv.args)
        
        # Model accuracy
        as_nuevo = mean(sapply(cv.mse_new, '[[', 'loss'))
        
        
        # CHECK IF THE NEW MODEL OUTPERFORMS THE PREVIOUS ONE
        if(maximize){
          improve = as_nuevo > as_actual
        }else{
          improve = as_nuevo < as_actual
        }
        if(improve){
          bestSet = varsModeloNew
          outvar = colnames(data)[varsModelo[j]]
          
          as_actual = as_nuevo
          cv.mse = cv.mse_new
        }
      }
      if(verbose){
        cat("ADDITION \n")
      }
    }

    # ADDITION
    ###############-
    varsModeloNew = c(varsModelo,newVar)
    dfModelo = data[,c(varsModeloNew, target_index)]

    
    # FIT NEW MODEL
    motbf.cv.args = append(list(data = dfModelo), args)
    cv.mse_new = do.call(motbf.cv, motbf.cv.args)
    
    # Model accuracy
    as_nuevo = mean(sapply(cv.mse_new, '[[', 'loss'))
    
    # CHECK IF THE NEW MODEL OUTPERFORMS THE PREVIOUS ONE
    if(maximize){
      improve = as_nuevo > as_actual
    }else{
      improve = as_nuevo < as_actual
    }
    if(improve){
      outvar = NULL
      bestSet = varsModeloNew
  
      as_actual = as_nuevo
      cv.mse = cv.mse_new
    }
    
    
    if(!is.null(bestSet)){
      varsModelo = bestSet
      globalImprove = TRUE
      
      if(verbose){
        if(!is.null(outvar)){ # replacement has improved more that addition
          
          cat("* Updated model:", colnames(data)[newVar], "REPLACES", outvar,"\n",
              " Variables included in the current model: ", paste0(colnames(data[,c(varsModelo, target_index)]), collapse = ', '),
              "\n ", loss_info ,"=", as_actual, "\n\n")
        }else{
          
          cat("* Updated model, variables:", paste0(colnames(data[,c(varsModelo, target_index)]), collapse = ', '), "\n ", loss_info ,"=", as_actual, "\n\n")
          
        }
      }
    }else{
      if(verbose){
        cat("- Reject variable:", colnames(data)[newVar], "\n ", loss_info ,"=", as_nuevo, "\n\n")
      }
    }
    
    
  }# end for loop over ordered set of explanatory variables
    if(method %in% c('gf', 'iwsr')){
      if(verbose){
        cat("Improve:", globalImprove,"\n" )
      }
    }else if(method == 'forward'){
      globalImprove = FALSE
    }

  }# end while loop
  
  if(verbose){
    cat("** Best model: ", paste0(colnames(data[,c(varsModelo, target_index)]), collapse = ', '), ".\n** ",loss_info, ": ", as_actual, "\n\n", sep = '')
  }
  dfModelo <- data[, sort(c(varsModelo, target_index))]
  
  return(list(data = dfModelo, index = c(varsModelo, target_index), loss = as_actual, crossvalidation = cv.mse))

}



