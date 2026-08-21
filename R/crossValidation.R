#' Cross-validation for MoTBFs
#' 
#' Perform a k-fold cross validation for Bayesian networks of class MoTBF.
#' 
#' If the basis function is MTE, the loss function is the log-likelihood of the model. 
#' On the other hand, if the basis function if MOP, the loss function might be either 
#' the predictive error ('pred') or the the log-likelihood of the model ('logl').
#' @param data  an object of class \code{"data.frame"}, which can contain continuous and discrete variables.
#' @param dag a network of the class \code{"bn"}, \code{"graphNEL"} or \code{"network"}.
#' @param k an integer indicating the number of folds to split the data set. If k =  0, the train and test sets
#' are the same data; if k = 1, hold-out validation is carried out, i.e., the data set is split in train (80% by default) and test (20% by default);
#' finally, if k >=2, k-fold cross validation is carried out.
#' @param loss a character string indicating which loss function should be used. 
#' Currently, two options are available: 'logl', for the log-likelihood of the model; and 'pred', for the predictive error. See details.
#' @param target a character string indicating which node is the target. This argument might be NULL if
#' the loss function chosen is the log-likelihood of the model ('logl').
#' @param fit.args a list containing optional arguments used to fit the models. These arguments must be those accepted by function 
#' \link{motbf.fit}, i.e., 'numIntervals' (4), 'POTENTIAL_TYPE' ('MOP'), 'maxParam' (NULL), 's' (NULL), 'priorData' (NULL) or 'scale' (TRUE). 
#' If fit.args is left NULL, the default values (in brackets) for those arguments will be used.
#' @param seed an integer to specify the seed. 
#' The k-folds are created randomly, so one might expect slightly different results unless 'seed' is used.
#' @param loss.args a list containing optional arguments related to the loss functions. Currently available arguments are:
#'  'loss.matrix' and 'percentage_test'. See details.
#' @param foldsIndex a list containing the test indexes for each fold of cross-validation. If is not null, \code{k} and \code{seed} are ignored.
#' @param ... Additional arguments. Not used currently.
#' @details 
#' Details on the \code{loss} argument:
#' \describe{
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
#' Note that the diagonal of the loss matrix is regarded as a reward, while the off-diagonal is regarded as a cost.
#' }
#' }
#' @return An object of class \code{"motbf.fit.cv"}. This is a list of \code{"k"} elements, 
#' each of them containing the results of each fold. More specifically, each element contains:
#' \describe{
#'  \item{test}{a \code{"data.frame"} of the subset used to fit the model.}
#'  \item{fitted}{an object of class \code{"motbf.fit"}, i.e., the model fitted.}
#'  \item{loss}{the loss value computed for the fold.
#'    If \code{"pred"} is chosen for the argument \code{"loss"}, then each element of the \code{"motbf.fit.cv"} object
#'    also contains:
#'    \describe{
#'      \item{predicted}{a vector of predictions for the target variable.}
#'      \item{observed}{a vector of observed records of the target variable.}
#'      }
#'    }
#' }
#' @examples 
#' #################
#' ### EXAMPLE 1 ###
#' #################
#' ## Perform 2-fold cross validation using the default model arguments
#' ## and the log-likelihood as loss function
#' 
#' # Load data
#' data(ecoli)
#' ecoli <- ecoli[,-c(1,9)]
#' 
#' # Learn DAG
#' dag <- LearningHC(ecoli)
#' 
#' # Run cross validation
#' cv = motbf.cv(data = ecoli, dag, k = 2, loss = 'logl')
#' cv
#' 
#' \donttest{
#' #################
#' ### EXAMPLE 2 ###
#' #################
#' ## Choose different arguments to fit the model parameters
#' fit.args = list(numIntervals = 3, POTENTIAL_TYPE = 'MOP', maxParam = 4)
#' 
#' # Run cross validation using the classification accuracy as loss function
#' cv = motbf.cv(data = ecoli, dag, k = 2, loss = 'pred', target = 'lip', 
#'   fit.args = fit.args)
#' cv
#' summary(cv)
#' 
#' #################
#' ### EXAMPLE 3 ###
#' #################
#' 
#' ## Specify a loss matrix to increase the penalty of classification errors
#' lossFunctionMatrix = matrix(c(c(1,2), c(3, 1)),nrow = 2, ncol = 2, byrow = TRUE)
#' 
#' # Run cross validation using the weighted classification accuracy as loss function
#' cv = motbf.cv(data = ecoli, dag, k = 2, loss = 'pred', target = 'lip', 
#'   loss.args = list(loss.matrix = lossFunctionMatrix))
#' cv
#' summary(cv)
#'   
#' #################
#' ### EXAMPLE 4 ###
#' #################
#' 
#' # Run cross validation using the root mean squared error as loss function
#' cv = motbf.cv(data = ecoli, dag, k = 2, loss = 'pred', target = 'mcg')
#' cv
#' }
#' @export
motbf.cv <- function(data, dag, loss, k = 10, target = NULL, seed = NULL, 
                     fit.args = NULL, loss.args = NULL, foldsIndex = NULL, ...){
  # browser()

  # if(!(loss %in% c('logl', 'rmse', 'accuracy', 'recall', 'precision', 'fscore', 'bscore', 'gmean'))){
  #   stop("Argument 'loss' must be one of the following: 'logL', 'rmse', 'accuracy', 'recall', 'precision', 'fscore', 'bscore', 'gmean'.")
  # }
  
  if(!(loss %in% c('logl', 'pred'))){
    stop("Argument 'loss' must be one of the following: 'logl', 'pred'.")
  }
  if(!is.null(fit.args$POTENTIAL_TYPE) && fit.args$POTENTIAL_TYPE == 'MTE' & loss == 'pred'){
    stop("The computation of the predictive error is not available for potentials of class MTE. Use loss = 'logl' or POTENTIAL_TYPE == 'MOP'.")
  }
  
  for (i in 1:length(data)) {
    if(is.character(data[,i])){
      data[,i]=as.factor(data[,i])
    }
  }
  
  if(loss == 'pred' & is.null(target)){
    stop("Use the argument 'target' to provide the variable of interest.")
  }

  if(!is.null(target) & !any(colnames(data)%in%target)){
    stop('Target name given is not in the dataset. Check name.')
  }
  
  if(loss == 'logl'){# COMPUTE MODEL LOGLIKELIHOOD
    lossfunction = logLikelihood.MoTBFBN
    lossfunction.args = list(object = NULL, data = NULL)
    
    target.type = NULL
  }else{# PREDICT TARGET VALUE
    lossfunction = stats::predict
    lossfunction.args = list(object = NULL, target = target, data = NULL, prob = TRUE)
    
    target.numeric = is.numeric(data[[target]])
    target.type = ifelse(target.numeric, 'Continuous', 'Discrete')
  }
  # Compute k-folds
  if(is.null(foldsIndex)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    # Split dataset in k-train and test sets
    if(k == 0){
      folds = list(list(Training = data, Test = data))
      K = 1
    }else if(k == 1){
      if(is.null(loss.args$p.test)){
        p.test = 0.2
      }else{
        p.test = loss.args$p.test
      }
      
      folds = list(TrainingandTestData(data, percentage_test = p.test))
      K = 1
    }else{
      folds = splitFolds(data,k)
      K = k
    }
  }else{# Se han introducido los indices de test de cada fold
    K = length(foldsIndex)
    folds <- list()
    for(i in 1:K){
      folds[[i]] <- list(Training = data[-foldsIndex[[i]],],
                         Test = data[foldsIndex[[i]],])
    }
  }
  
  
  # Save original folds (object folds might be scaled, or missing observations removed)
  fold_original = folds
  scale = FALSE
  if(is.null(fit.args$scale)||fit.args$scale == TRUE){
    scale = TRUE # save for later
    # save vector of means and sd to scale the training and test sets in each fold
    v_mean = sapply(data, function(x){ifelse(is.numeric(x), mean(x, na.rm = TRUE),NA)})
    v_sd = sapply(data, function(x){ifelse(is.numeric(x), sd(x, na.rm = TRUE),NA)})
    
    for(i in 1:length(folds)){
      folds[[i]]$Training = rescale_data(folds[[i]]$Training, v_mean = v_mean, v_sd = v_sd)
      folds[[i]]$Test = rescale_data(folds[[i]]$Test, v_mean = v_mean, v_sd = v_sd)
    }

    if(is.null(fit.args)){
      fit.args = list(scale = FALSE)
    }else{
      fit.args$scale = FALSE
    }
  }
  
  results = list()
  res = list(train = NULL, test = NULL, fitted = NULL, loss = NULL)
  # Perform k-fold cross validation
  for(i in 1:K){
    # Usar los datos train y test para aprender (train) y predecir (test) los modelos 
    trainData <- folds[[i]]$Training # trainData is checked inside motbf.fit()
    testData <- folds[[i]]$Test
    testData = check_data(testData)# missing observations might be removed
    
    
    if(is.character(dag)){
      # Get model structure
      if(!(dag %in% c('NB', 'TAN', 'HC'))){
        stop("Argument 'dag' must be an object of class 'bn' or a character string: 'NB' (naive Bayes structure); 'TAN' (tree augmented Naive Bayes structure); or 'HC' (structure via hill-climbing algorithm).")
      }
      dag = getStructure(trainData, method = dag, target = target)
    }
    # APRENDER los parametros de la red
    motbf.fit.args = append(list(graph = dag, data = trainData), fit.args)
    
    bn <- do.call(motbf.fit, motbf.fit.args)
    
    # PREDECIR
    lossfunction.args$object = bn
    lossfunction.args$data = testData
    
    pred = do.call(lossfunction, lossfunction.args)
    
    # SAVE results
    res$train <- fold_original[[i]]$Training
    res$test <- fold_original[[i]]$Test
    
    if(scale){
      # res$train <- fold_original[[i]]$Training
      # res$test <- fold_original[[i]]$Test
      
      potential = ifelse(is.null(fit.args$POTENTIAL_TYPE), 'MOP', fit.args$POTENTIAL_TYPE)
      # Return fitted model in original scale
      res$fitted = rescale_motbf_fit(bn, POTENTIAL_TYPE = potential, data = data )
    }else{
      # res$train <- trainData
      # res$test <- testData
      res$fitted <- bn
    }
    
    if(loss == 'logl'){
      res$loss = pred
    }else{
      obs = testData[[target]]
      
      if(target.numeric){# If target variable is numeric, compute the RMSE
        if(scale){
          m = mean(data[[target]])
          s = sd(data[[target]])
          pred = pred*s+m
          obs = obs*s+m
        }
        loss.value = rmse(obs, pred)
      }else{# If target variable is factor, compute the classification error
        loss.matrix = loss.args[['loss.matrix']]
        loss.value = accuracy(obs, pred, loss.matrix)
      }
      res$loss = loss.value
      res$predicted = pred
      res$observed = obs
    }

    attr(res$test, 'index') = attr(folds[[i]], 'index')$test
    attr(res$train, 'index') = attr(folds[[i]], 'index')$train
    results[[i]] <- res

  }
  # browser()
  attr(results, 'loss.info') = loss
  attr(results, 'k') = k
  attr(results, 'target') = target
  attr(results, 'target.type') = target.type
  attr(results, 'loss.matrix') = loss.args['loss.matrix']
  results = new_motbf_fit_cv(results)
  
  return(results)
}

#' Confusion Matrix
#'
#' Computes the confusion matrix for a given fitted model object.
#'
#' @param x An object used to select the appropriate method, such as one of class \code{'motbf_fit_cv'}.
#' @param digits An integer indicating the number of decimal places to round the results.
#' @param ... Further arguments. Not used currently.
#'
#' @return A matrix or table containing the confusion matrix.
#' @export
confusionMatrix <- function(x, ...){
  UseMethod('confusionMatrix')
}
#' @rdname confusionMatrix
#' @exportS3Method
confusionMatrix.motbf_fit_cv <- function(x, digits = 4,...){
 if(!is.motbf_fit_cv(x)){
   stop("Object 'x' must be of class 'motbf_fit_cv'. See function motbf.cv()")
 }
  target.type = attr(x, 'target.type')
  
  if(is.null(target.type) | target.type != 'Discrete'){
    stop("Target variable must be of type discrete.")
  }
  k = length(x)
  
  loss.matrix = attr(x, 'loss.matrix')[[1]]
  
  # Crear vectores para almacenar resultados
  ns <- ncol(x)
  fold_res <- list()
  # cat("Observed values on rows\nPredicted values on columns\n")
  for(i in 1:k){
    if(k >1){
      cat("\n--------- Fold", i, "---------\n")
    }else if(attr(x, 'k') == 1){
      cat("\n--------- Hold-out validation ---------\n")
    }else{
      cat("\n--------- Validation using training set ---------\n")
    }
    
    z = x[[i]]
    cm = table(z$observed, z$predicted)
    names(dimnames(cm)) <- c('Observed', 'Predicted')
    print(cm)
    cat("\n")
    st = colnames(cm)
    
    accuracy <- sum(diag(cm))/sum(cm)
    
    precision <- diag(cm) / colSums(cm)# TP / (TP+FP)
    names(precision) = st
    
    recall <- (diag(cm) / rowSums(cm))# TP / (TP+FN)
    names(recall) = st
    
    fscore <- (2*precision*recall)/(precision+recall)
    names(fscore) = st
    
    gmean = gm_mean(recall)
    

    fold_res[[i]] = list(accuracy = accuracy, precision = precision, recall = recall, fscore = fscore, gmean.recall = gmean)
    print(t(as.data.frame( fold_res[[i]][2:4])))
    cat("Accuracy =", round(accuracy, digits), "\n")
    cat("Geometric mean of recalls = ", round(gmean, digits), '\n')
    prob = attr(z$predicted, 'prob')
    if(!is.null(prob)){
        comp.BS <- as.data.frame(prob)
        comp.BS$OBS <- z$observed

        BS <- brier.score(comp.BS)
        fold_res[[i]] = append(fold_res[[i]], list(BS = BS))
        cat("Brier score =", round(BS, digits), "\n")
    }
    if(!is.null(loss.matrix)) {
      # weighted accuracy
      wAccuracy = sum(diag(as.matrix(cm)) * diag(loss.matrix))/sum(cm*loss.matrix) 
      # sum(w*m)/sum(m); w= loss.matrix; m = cm
      fold_res[[i]] = append(fold_res[[i]], list(weighted.accuracy = wAccuracy))
      cat("Weighted accuracy =", round(wAccuracy, digits), "\n")
    }
  }# end loop over folds
  if(k >1){
    cat("\n-----------------------------------")
    cat("\n--------- Average metrics ---------\n")
    
    
    av_accuracy = round(mean(sapply(fold_res, '[[', 'accuracy')), digits)
    av_precision = round(Reduce('+',lapply(fold_res, '[[', 'precision'))/k, digits)
    av_recall = round(Reduce('+',lapply(fold_res, '[[', 'recall'))/k, digits)
    av_fscore = round(Reduce('+',lapply(fold_res, '[[', 'fscore'))/k, digits)
    
    av_gmean= gm_mean(av_recall)
    
    global_res = list(accuracy = av_accuracy, precision = av_precision, recall = av_recall, fscore = av_fscore, gmean.recall = av_gmean)
    
    print(t(as.data.frame( global_res[2:4])))
    cat("Accuracy =", round(av_accuracy, digits), "\n")
    cat("Geometric mean of recalls = ", round(av_gmean, digits), '\n')
    
    if(!is.null(prob)){
      av_BS = round(mean(sapply(fold_res, '[[', 'BS')), digits)
      global_res = append(global_res, list(BS = av_BS))
      cat("Brier score =", round(av_BS, digits), "\n")
    }
    if(!is.null(loss.matrix)) {
      # weighted accuracy
      wAccuracy = sum(diag(as.matrix(cm)) * diag(loss.matrix))/sum(cm*loss.matrix) 
      av_wAccuracy = round(mean(sapply(fold_res, '[[', 'weighted.accuracy')), digits)
      global_res = append(global_res, list(weighted.accuracy = av_wAccuracy))
      cat("Weighted accuracy =", round(av_wAccuracy, digits), "\n")
    }
    results = list(global = global_res, fold = fold_res)
  }else{
    results = list(global = fold_res)
  }
  
return(invisible(results))
}


#' @importFrom stats cor median qqplot
#' @importFrom graphics hist boxplot lines
#' @noRd
gof <- function(x, digits = 4, method, plot = FALSE){
  if(!is.motbf_fit_cv(x)){
    stop("Object 'x' must be of class 'motbf_fit_cv'. See function motbf.cv()")
  }
  target.type = attr(x, 'target.type')
  target = attr(x, 'target')
  
  if(method == 'pred' & (is.null(target.type) || target.type != 'Continuous')){
    stop("Target variable must be of type continuous")
  }
  
  
  
  opar <- par(no.readonly =TRUE)       
  on.exit(par(opar)) 
  
  k = length(x)
  
  par(mfrow = c(2,2))
  fold_res <- list()
  df = data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c('Observed', 'Predicted', 'Fold'))))
  for(i in 1:k){
    if(k >1){
      cat("\n--------- Fold", i, "---------\n")
    }else if(attr(x, 'k') == 1){
      cat("\n--------- Hold-out validation ---------\n")
    }else{
      cat("\n--------- Validation using training set ---------\n")
    }
    
    z = x[[i]]
    if(method == 'logl'){
      logl = z$loss
      fold_res[[i]] = list(logl = logl)
      row_names = 'Log-likelihood'
      target = ' '
    }else{
      df_i = data.frame(Observed = z$observed, Predicted = z$predicted, Fold = i)
      df = rbind(df, df_i)
      
      
      y = z$observed
      y_fit = z$predicted
      
      rmse = z$loss
      # Normalized root mean squared error
      n_mean = rmse/mean(y)
      # n_sd = rmse/sd(y)
      # n_maxmin = rmse/diff(range(y))
      # n_iqr = rmse/IQR(y)
      
      correlation = cor(y, y_fit)
      # browser()
      # R-squared
      r.squared = 1-sum((y-y_fit)^2)/sum((y-mean(y))^2)
      # Adjusted R-squared
      n = length(y) # number of samples
      p = ncol(z$train)-1 # number of predictors
      adj.r.squared = 1-(1-r.squared)*(n-1)/(n-p-1)
      
      # residuals
      e = y_fit - y
      # Bias
      bias = mean(e)
      # Mean absolute error
      mae = mean(abs(e))
      # Normalized mean absolute error
      n_mae = mae/mean(y)
      # Mean absolute percentage error
      mape = mean(abs(e)/y)*100
      # Median absolute deviation
      mad = median(abs(e-median(e)))
      
      # naive forecast
      d_j1 = z$train[[target]]
      d_j = d_j1[-1]
      N = length(d_j1)
      d_j1 = d_j1[-N]
        # Q is the naive forecast computed on the training data
        # MAE is calculated on the test data and Q is calculated on the training data
      Q = sum(abs(d_j-d_j1))/(N-1)
      # Mean Absolute Scaled Error
      mase = mae/Q
        # If MASE is less than 1, it means that the forecast is better than the one step naive method on training data. 
        # If it is more than 1, it means that the forecast method is worse than the forecast using one step naive forecasting approach on the training data.
      
      
      # fold_res[[i]] = list(rmse = rmse, nrmse_sd = n_sd, nrmse_mean = n_mean, nrmse_maxmin = n_maxmin, nrmse_iqr = n_iqr, correlation = correlation)
      fold_res[[i]] = list(rmse = rmse, nrmse = n_mean, correlation = correlation, r.squared = r.squared, adj.r.squared = adj.r.squared, mae = mae, nmae = n_mae, mape = mape, mad = mad, mase = mase)
      
      row_names = c('Root mean squared error (RMSE)', 
                    'Normalized RMSE (by mean)',
                    'Correlation',
                    'R-squared',
                    'Adjusted R-squared',
                    'Mean absolute error (MAE)',
                    'Normalized MAE (by mean)',
                    'Mean absolute percentage error (MAPE)',
                    'Median absolute deviation (MAD)',
                    'Mean Absolute Scaled Error (MASE)')
    }
    
    a = t(as.data.frame(fold_res[[i]]))
    colnames(a) = target
    rownames(a) = row_names
    print(a, digits = digits)
    
}
  if(k >1){
    # browser()
    
    
    cat("\n-----------------------------------")
    cat("\n--------- Average metrics ---------\n")
    
    if(method == 'logl'){
      av_logl = round(mean(sapply(fold_res, '[[', 'logl')), digits)
      global_res = list(logl = av_logl)
      
    }else{
      global_res = as.list(round(Reduce('+', lapply(fold_res, unlist))/k, digits))
      # av_rmse = round(mean(sapply(fold_res, '[[', 'rmse')), digits)
      # av_n_sd = round(mean(sapply(fold_res, '[[', 'nrmse_sd')), digits)
      # av_n_mean = round(mean(sapply(fold_res, '[[', 'nrmse_mean')), digits)
      # av_n_maxmin = round(mean(sapply(fold_res, '[[', 'nrmse_maxmin')), digits)
      # av_n_iqr = round(mean(sapply(fold_res, '[[', 'nrmse_iqr')), digits)
      # av_corr = round(mean(sapply(fold_res, '[[', 'correlation')), digits)
      # 
      # global_res = list(rmse = av_rmse, nrmse_sd = av_n_sd, nrmse_mean = av_n_mean, nrmse_maxmin = av_n_maxmin, nrmse_iqr = av_n_iqr, correlation = av_corr)
    }
    
    a = t(as.data.frame(global_res))
    colnames(a) = target
    rownames(a) = row_names
    print(a, digits = digits)
    
    results = list(global = global_res, fold = fold_res)
  }else{
    
    global_res = as.list(unlist(fold_res))
    results = list(global = fold_res)
  }

  if(method == 'pred'){
    if(plot %in% "univ" | plot == TRUE){
      ## PLOT 1 ##
      plot(df$Observed, df$Predicted, pch = 20, xlab = 'Observed', ylab = 'Predicted', main = 'Observed vs. predicted values')
      
      ## PLOT 2 ##
      lim = c(floor(min(df$Observed, df$Predicted)), ceiling(max(df$Observed, df$Predicted)))
      qqplot(df$Observed, df$Predicted, pch = 20, xlab = 'Observed', ylab = 'Predicted', main = 'Q-Q plot Obs. vs. Pred.', 
           xlim = lim, ylim = lim)
      abline(a = 0, b = 1, col = '#ffc425')
      # plot(df$Observed-df$Predicted, pch = 19, ylab = 'Residuals', main = 'Model residuals')
      # abline(h = mean(df$Observed-df$Predicted), col = '#ffc425', lwd = 2)
      
      ## PLOT 3 ##
      # hist(df$Observed, main = 'Histogram of observed target')
      # hist(df$Predicted, main = 'Histogram of predicted target', add = TRUE)
      h = max(c(hist(df$Observed, plot = FALSE)$counts, hist(df$Predicted, plot = FALSE)$counts))
      hist(df$Observed, prob = FALSE, col = "#FFC107", density = 10, angle = -45, xlab = target, main = paste('Histogram of', target), ylim = c(0, h*1.5))
      hist(df$Predicted, prob = FALSE, col = "#0C7BDC", main = "",  cex.lab=1.5, cex.axis=1.5, density = 10,add = TRUE)
      legend('topleft', legend = c("Observed", "Predicted"), fill = c("#FFC107", "#0C7BDC"), density = 20, angle = c(-45, 45), horiz = TRUE)
      
      ## PLOT 4 ##
      boxplot(df[,c(1,2)])
      # val = unlist(fold_res)
      # index_nrmse = grep('_', names(val))
      # index = (1:(6*k))[!(1:(6*k)%in%index_nrmse)]
      # 
      # ## PLOT 3 ##
      # ## unnormalized rmse and correlation
      # col = rep(gg_color_hue(2), times = k)
      # plot(rep(1:k, each = 2), val[index], col = col, ylim = c(0, max(val[index])*1.5), pch = 20,
      #      main = "Unnormalized RMSE and correlation",
      #      xlab = 'Fold', ylab = 'Value')
      # 
      # legend("topright", legend = c("rmse" , "cor"), 
      #        col = unique(col), cex = 0.8, lty = 1, ncol = 2)
      # abline(h = global_res[c(1,6)], col = col, lty = 2)
      # 
      # ## PLOT 4 ##
      # ## normalized rmse 
      # col = rep(gg_color_hue(4), times = k)
      # plot(rep(1:k, each = 4), val[index_nrmse], col = col, ylim = c(0, max(val[index_nrmse])*1.5), pch = 20,
      #      main = "Normalized RMSE",
      #      xlab = 'Fold', ylab = 'Value')
      # 
      # legend("topright", legend = c("nrmse_sd", "nrmse_mean", "nrmse_maxmin", "nrmse_iqr"), 
      #        col = unique(col), cex = 0.8, lty = 1, ncol = 2)
      # abline(h = global_res[c(2:5)], col = col, lty = 2) 
      
    } 
    if(plot %in% "ts" | plot == TRUE){
      ## Plot the temporal data series
      if(attr(x, 'k') == 0){
        par(mfrow = c(1,1))
        plot(df$Observed, type = 'l', col = 'red', xlab = 'Time', ylab = 'Value')
        lines(df$Predicted, type = 'l', col = 'blue')
        legend('topleft', legend = c('Observed', 'Predicted'),
               col = c('red', 'blue'), lwd = 2)
      }else{
        index = unlist(lapply(lapply(x, '[[', 'test'), attr, 'index'))
        
        target = attr(x, 'target')
        
        pred = unlist(lapply(x, '[[', 'predicted'))
        obs = unlist(lapply(lapply(x, '[[', 'test'), '[[', target))
        
        
        par(mfrow = c(1,1))
        plot(obs[order(rank(index))], type = 'l', col = 'red', xlab = 'Time', ylab = 'Value')
        lines(pred[order(rank(index))], col = 'blue')
        legend('topleft', legend = c('Observed', 'Predicted'), 
               col = c('red', 'blue'), lwd = 2)
      }
    }
 
    


  }
  par(opar)
  return(results)
}

gm_mean = function(a){prod(a)^(1/length(a))}



#Función para calcular el Brier Score. El argumento de entrada es un dataframe
#donde las columnas son las probabilidades con las que se predicen cada uno de los estados de 
#la variabe objetivo, siendo la última de éstas el verdadero valor de la variable.
brier.score <- function(df){
  
  s <- ncol(df)
  # Observed outcome
  obs = df[,s]
  # Predicted probability distribution
  pd = df[,-c(s)]
  
  # Check names on 'pd'. If state names are 'numbers' data.frame adds an 'X' at the begining.
  st = names(pd)
  if(any(!(st%in%levels(obs)))){
    st = gsub('X', '', st)
    if(all(!(st%in%levels(obs)))){
      stop("Names in probability distribution do not match variable levels.")
    }
  }
  
  N = length(obs)
  
  # Position of predicted class (max probability outcome)
  pos_st_max = apply(pd, 1, which.max)
  
  # Create one-hot encoding for observed state
  pos_out = sapply(obs, function(x){which(x==st)})
  
  # Initialize data.frame
  st_onehot = data.frame(matrix(0, ncol = length(st), nrow = N, dimnames = list(NULL, st)))
  for(i in 1:N){
    st_onehot[i,pos_out[i]] = 1
  }
  
  # Compute Brier Score: 1/N*sum_i^N sum_k^st (pd_ik-output_ik)^2; where output_ik is coded as 0 and 1
  BS = sum((pd-st_onehot)^2)/N
  return(BS)
  
  # Old code
  # BS <- c()
  # for(j in 1:nrow(df)){
  #   
  #   Obs <- rep(0,s-1)
  #   mp<-which.max(df[j,1:(s-1)])
  #   
  #   if(st[mp]==df[j,s]){
  #     Obs[mp]=1
  #   }
  #   
  #   BS[j] = sum((df[j,1:(s-1)] - Obs)^2)
  # }
  # return(mean(BS))
}


rmse <- function(x, y){
  sqrt(mean((x-y)^2, na.rm = TRUE))
}

accuracy <- function(x, y, loss.matrix = NULL){
  # x: observed values
  # y: predicted values
  # browser()
  cm <- table(x, y)
  
  # Calcula las metricas para evaluar el clasificador
  if(is.null(loss.matrix)){
    accuracy <- sum(diag(cm))/sum(cm)
  }else{
    accuracy = sum(diag(as.matrix(cm)) * diag(loss.matrix))/sum(cm*loss.matrix)
  }
  
  
  
  
  return(accuracy)
}