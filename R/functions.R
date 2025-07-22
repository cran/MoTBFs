#' Data pre-processing utilities
#' 
#' Collection of functions for discretizing, standardizing, converting factors to 
#' characters and other usufull methods for pre-processing datasets.
#' 
#' @name dataMining
#' @rdname dataMining
#' @param dataset A dataset of class \code{"data.frame"}. Tha variables of the dataset can be discrete and continuous.
#' @param discreteVariables A \code{"character"} array with the names of the discrete variables.
#' @param numIntervals Number of bins used to discretize the continuous variables.
#' @param factor A boolean value indicating if the variables should be considered as
#' \code{"factor"} or as \code{"character"}. By default it is set to \code{FALSE}.
#' @param binary By default it is set to \code{FALSE}, indicating that only binary entries are 
#' used for continuous variables; a \code{TRUE} value means that binary entries are used to 
#' discretize the full dataset taking into account the states the discrete variables.
#' @param discreteData A discretized dataset of class \code{"data.frame"}.
#' @param DiscreteVariablesStates The output of the function \code{discreteVariablesStates}.
#' @param namevariables an array with the names of the varibles.
#' @param scale A \code{"numeric"} vector (when it refers to a single variable) or a \code{"list"} 
#' containing the name(s) of the variable(s) and the scale value.
#' @param X A \code{"numeric"} vector with the data values of a continuous variable.
#' @details
#' \code{whichDiscrete()} selects the position of the discrete variables.
#' 
#' \code{discreteVariables_as.character()} transforms the values of the discrete variables into character values.
#' 
#' \code{standardizeDataset()} standardizes all the variables in a data set.
#' 
#' \code{discretizeVariablesEWdis()} discretizes the continuous variables in a dataset using 
#' equal width binning.
#' 
#' \code{discreteVariablesStates()} extracts the states of the qualitative variables.
#' 
#' \code{nstates()} computes the number of different values of the discrete variables.
#' 
#' \code{quantileIntervals()} gets the quantiles of a variable taking into account the number of intervals
#' into which its domain is splitted.
#' @examples
#' ## dataset: 2 continuous variables, 1 discrete variable.
#' data <- data.frame(X = rnorm(100),Y = rexp(100,1/2), Z = as.factor(rep(c("s","a"), 50)))
#' disVar <- "Z" ## Discrete variable
#' class(data[,disVar]) ## factor
#' 
#' data <- discreteVariables_as.character(dataset = data, discreteVariables = disVar)
#' class(data[,disVar]) ## character
#' 
#' whichDiscrete(dataset = data, discreteVariables = "Z")
#' 
#' standData <- standardizeDataset(dataset = data)
#' 
#' disData <- discretizeVariablesEWdis(dataset = data, numIntervals = 3)
#' 
#' l <- discreteVariablesStates(namevariables = names(data), discreteData = disData)
#' 
#' nstates(DiscreteVariablesStates = l)
#' 
#' ## Continuous variables
#' quantileIntervals(X = data[,1], numIntervals = 4)
#' quantileIntervals(X = data[,2], numIntervals = 10)
#' @export
whichDiscrete <- function(dataset, discreteVariables) which(colnames(dataset)%in%discreteVariables)

#' @rdname dataMining
#' @export
discreteVariables_as.character <- function(dataset, discreteVariables)
{
  if(is.numeric(discreteVariables)) pos <- discreteVariables
  else pos  <- whichDiscrete(dataset,discreteVariables)
  for(i in 1:ncol(dataset)){
    if(any(pos==i)) dataset[,i]  <- as.character(dataset[,i])
    else dataset[,i] <- as.numeric(as.character(dataset[,i]))
  }
  return(dataset)
}

#' @rdname dataMining
#' @export
standardizeDataset <- function(dataset)
{
  for(i in 1:length(dataset)) if(is.numeric(dataset[,i])) dataset[,i] <- scale(dataset[,i])[,1]
  return(as.data.frame(dataset))
}

#' @rdname dataMining
#' @export
discretizeVariablesEWdis <- function(dataset, numIntervals, factor=FALSE, binary=FALSE)
{
  for(i in 1:ncol(dataset)) {
    if(is.character(dataset[,i])){
      if(!binary) dataset[,i] <- dataset[,i]
      if(binary){
        states <- discreteVariablesStates(colnames(dataset)[i],dataset)
        nint <- nstates(states)
        for(j in 1:nint) dataset[dataset[,i]==states[[1]]$states[j],i] <- j-1
      }
    }
    else dataset[,i] <- as.character(as.vector(cut(as.numeric(dataset[,i]), numIntervals,labels=as.numeric(0:(numIntervals-1)))))
  }
  if(factor) for(i in 1:ncol(dataset)) dataset[,i] <- as.factor(dataset[,i])
  
  return(dataset)
}

#' @rdname dataMining
#' @export
discreteVariablesStates <- function(namevariables, discreteData)
{
  variablesStates <- c()
  for(i in 1:length(namevariables)){
    states <- sort(discreteData[,namevariables[i]][which(duplicated(discreteData[,namevariables[i]])==F)])
    info <- list(variable=namevariables[i], states=states)
    variablesStates[[length(variablesStates)+1]] <- info
  }
  return(variablesStates)
}

#' @rdname dataMining
#' @export
nstates <- function(DiscreteVariablesStates) sapply(1:length(DiscreteVariablesStates), 
                                                    function(i) length(DiscreteVariablesStates[[i]]$states))

#' @rdname dataMining
#' @export
quantileIntervals <- function(X, numIntervals) sapply(1:(numIntervals-1), function(i) quantile(X, i/numIntervals))

#' @rdname dataMining
#' @export
scaleData <- function(dataset, scale)
{
  if(is.numeric(dataset)) return(dataset*scale[[1]])
  else{
    names <- names(scale)
    for(i in 1:length(names)){
      dataset[,names[i]] <- dataset[,names[i]]*scale[[i]]  
    }
    return(dataset)
  }   
}

#' Dataset subsetting
#'
#' Collection of functions for subsetting a \code{"data.frame"} by rows or columns, and
#' to create training and test partitions.
#' 
#' @name subsetData
#' @rdname subsetData
#' @param data A dataset of class \code{data.frame}.
#' @param percentage_test The proportion of data that goes to the test set (between 0 and 1).
#' @param discreteVariables A \code{character} vector with the name of the discrete variables.
#' @param nameX A \code{character} vector with the name of the child variable in the conditional method.
#' @param nameY A \code{character} vector with the name of the parent variables in the conditional method.
#' @param nameVariable A \code{character} vector with the name of the variable to be filtered.
#' @param min,max Boundary values to filter out.
#' @return \code{TrainingandTestData()} returns a list of 2 elements containing the train and test datasets. 
#' \code{newData()} and \code{splitdata()} return a subset of variables or observations, respectively.
#' @examples
#' \donttest{
#' ## Dataset
#' X <- rnorm(1000)
#' Y <- rchisq(1000, df = 8)
#' Z <- rep(letters[1:10], times = 1000/10)
#' data <- data.frame(X = X, Y = Y, Z = Z)
#' data <- discreteVariables_as.character(dataset = data, discreteVariables ="Z")
#' 
#' ## Training and Test Datasets
#' TT <- TrainingandTestData(data, percentage_test = 0.2)
#' TT$Training
#' TT$Test
#' 
#' ## Subset Dataset
#' newData(data, nameX = "X", nameY = "Z")
#' splitdata(data, nameVariable = "X", min = 2, max= 3)
#' }
#' @export
TrainingandTestData <- function(data, percentage_test, discreteVariables=NULL)
{
  N <- nrow(data); n <- round(percentage_test*N)
  randomN <- sample(1:N,n)
  XTraining <- data[-randomN,]
  XTest <- data[randomN,]
  
  if(is.null(discreteVariables)) discreteVariables <- colnames(data)[sapply(data, is.character)]
  for(i in which(colnames(data)%in%discreteVariables)){
    trainingStates <- discreteVariablesStates(colnames(data)[i], XTraining)[[1]]$states
    states <- discreteVariablesStates(colnames(data)[i], data)[[1]]$states
    if(!all(states%in%trainingStates)){
      s <- states[!(states%in%trainingStates)]
      for(j in 1:length(s)){
        pos <- which(XTest[,i]==s[j])
        if(length(pos)>1) pos <- sample(pos,1)
        XTraining[nrow(XTraining)+1,] <- XTest[pos,]
        XTest <- XTest[-pos,]
      }
    }
  }
  
  for(i in which(!(colnames(data)%in%discreteVariables))){
    if(!(min(XTraining[,i])<=min(XTest[,i]))){
      pos <- which(min(XTest[,i])==XTest[,i])
      if(length(pos)>1) pos <- sample(pos, 1)
      XTraining[nrow(XTraining)+1,] <- XTest[pos,]
      XTest <- XTest[-pos,]
    }
    if(!(max(XTraining[,i])>=max(XTest[,i]))){
      pos=which(max(XTest[,i])==XTest[,i])
      if(length(pos)>1) pos <- sample(pos, 1)
      XTraining[nrow(XTraining)+1,] <- XTest[pos,]
      XTest <- XTest[-pos,]
    }
  }
  rownames(XTraining) <- 1:nrow(XTraining); rownames(XTest) <- 1:nrow(XTest)
  return(list(Training=XTraining, Test=XTest))
}

#' @rdname subsetData
#' @export
newData <- function(data, nameX, nameY) return(data[, c(nameX, nameY)])

#' @rdname subsetData
#' @export
splitdata <- function(data, nameVariable, min, max){
  parent <- data[, nameVariable]
  return(subset(data,((parent>=min)&(parent<=max))))
}

#' Data cleaning
#' 
#' Delete rows of a dataset wich contains anomalous values.
#' 
#' @param data A dataset of class \code{"matrix"} or \code{"data.frame"},
#' @param strangeElements A \code{"character"} string which contains the elementes to remove.
#' @export
preprocessedData <- function(data, strangeElements)
{
  for(j in 1:length(strangeElements)){
    n.row <- sapply(1:nrow(data), function(i) any(data[i,]==strangeElements[j]))
    data <- data[!n.row,]
  }
  return(data) 
}


#' Remove Objects from Memory
#' 
#' Clean the memory. Delete all the objects in memory and a garbage collection
#' takes place.
#' 
#' @param envir The currently active environment; by default It is the gloval environment.
#' @param n Number of garbage collection repetitions; by default \code{n = 2}.
#' @export
#' @examples
#' ## Run to clean the environment
#' clean()
#' clean(n=2)

clean <- function(envir = globalenv(), n = 2)
{
  ll <- ls(envir = envir)
  rm(list = ll, envir = envir)
  for(i in 1:n) gc()
}
