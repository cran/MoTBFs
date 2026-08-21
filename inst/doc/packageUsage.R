## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(MoTBFs)

## ----eval = FALSE-------------------------------------------------------------
# ## Install and load the MoTBFs package
#  install.packages("MoTBFs")
#  library("MoTBFs")

## -----------------------------------------------------------------------------
## Load ecoli dataset
 data("ecoli", package = "MoTBFs")

## -----------------------------------------------------------------------------
## Drop the first and last variables of the ecoli dataset 
  data <- ecoli[,-c(1,9)]

## ----eval = TRUE--------------------------------------------------------------
## Split the dataset into train and test subsets
  set.seed(2)
  dataTT <- TrainingandTestData(data, percentage_test = 0.2)
  trainingData <- dataTT$Training
  testData <- dataTT$Test

## -----------------------------------------------------------------------------
## Learn the structure of the Bayesian network using the training data
  dag <- LearningHC(trainingData)
  dag

## ----message=FALSE, eval = requireNamespace("Rgraphviz", quietly = TRUE)------
# Visualize the network structure in topological order
bnlearn::graphviz.plot(dag)

## -----------------------------------------------------------------------------
## Learn the density of variable mcg, using MTEs or MOP as basis functions
 f1 <- univMoTBF(trainingData[,1], POTENTIAL_TYPE = "MTE", nparam = 13)
 f2 <- univMoTBF(trainingData[,1], POTENTIAL_TYPE = "MOP", nparam = 11)

## -----------------------------------------------------------------------------
print(f2)

summary(f2)

## ----fig.align='center', fig.asp=1.1, out.width='80%'-------------------------
## Plot the densities f1 and f2 over the histogram of variable mcg
  hist(trainingData[,1], prob = TRUE , main = "", xlab = "X")
  plot(f1, xlim = range(trainingData[,1]), col = "red", add = TRUE, lwd = 3)
  plot(f2, xlim = range(trainingData[,1]), col = "blue", add = TRUE, lwd = 3, lty = 2)
  legend("topleft", legend = c("MTE", "MOP"),col = c("red", "blue"), lty = 1:2, lwd = 3, inset = c(0, -0.5), xpd = TRUE)

## -----------------------------------------------------------------------------
## Compute log-likelihood of the fitted densities
  sum(log(as.function(f1)(testData[,1])))

  sum(log(as.function(f2)(testData[,1])))

## Compute BIC score of the fitted densities  
  BICMoTBF(f1,testData[,1])

  BICMoTBF(f2,testData[,1])


## -----------------------------------------------------------------------------
## Simulate data from the estimated density f2
  set.seed(5)
  X <- rMoTBF(size = 400, fx = f2)
  
## Test whether or not the simulated sample and the observed data come from the same distribution
  ks.test(trainingData[,1], X)

## -----------------------------------------------------------------------------
## Compare the histogram of both distributions
 hist(trainingData[,1], prob = TRUE, col = "#FFC107", density = 10, angle = -45, main = "")
 hist(X, prob = TRUE, col = "#0C7BDC", main = "", ylim = c(0,2.2), cex.lab=1.5, cex.axis=1.5, density = 10,add = T)
 legend('topleft', legend = c("Training data", "Simulated data"), fill = c("#FFC107", "#0C7BDC"), density = 20, angle = c(-45, 45), inset = c(0, -0.5), xpd = TRUE)

## -----------------------------------------------------------------------------
## Compare the the CDF of both distributions
 plot(ecdf(trainingData[,1]), cex = 0, lwd = 3 , cex.lab = 1.5, cex.axis = 1.5, main = "")
 plot(integrate.motbf(f2), xlim = range(trainingData[,1]), col ="red", lwd = 3, add = TRUE)
 legend('topleft', legend = c("Training data", "Simulated data"), col = c("black", "red"), lwd = 3, inset = c(0, -0.5), xpd = TRUE)

## -----------------------------------------------------------------------------
## Compute the derivative and integral of the fitted density

  # Coefficients of the fitted density
  coef(f2)


  # Indefinite integral of the fitted density
  integrate.motbf(f2)

  
  # Definite integral of the fitted density
  integrate.motbf(f2, lower = min(trainingData[,1]), upper = max(trainingData[,1]))

  
  # Derivarive of the fited density
  derivMoTBF(f2)

## -----------------------------------------------------------------------------
## Learn joint distributions
  P = jointmotbf.fit(X = trainingData[,c("mcg", "alm1")], dimensions = c(5,5))
  
  attributes(P)


## -----------------------------------------------------------------------------
print(P)

summary(P)

## -----------------------------------------------------------------------------
## Plot the joint distribution of 2 variables
 par(mar=c(2,3,2,2))
 
 # Filled contour
 plot(P, data = trainingData[,c(1,6)])
 
 # Simple contour
 plot(P, data = trainingData[,c(1,6)], filled = FALSE, cex.lab = 2, cex.axis = 1.85, lwd = 1.5) 
  
 # Perspective
 plot(P, type = "perspective", data = trainingData[,c(1,6)], orientation=c(25,25), cex.lab = 2,xaxs = "i") 

## -----------------------------------------------------------------------------
marginal.jointmotbf(P, var = "mcg")
marginal.jointmotbf(P, var = 2)

## -----------------------------------------------------------------------------
## Learn conditional distributions
P <- conditionalMethod(trainingData, nameParents = "mcg", nameChild = "gvh",
                         numIntervals = 5, POTENTIAL_TYPE ="MOP", scale = FALSE)
printConditional(P)


## -----------------------------------------------------------------------------
par(mar=c(2,3,2,2))
## Plot the conditional density of gvh given mcg
plotConditional(P, data = trainingData, nameChild = "gvh", points = TRUE)

## ----message = FALSE----------------------------------------------------------
## Learn the distributions of a Bayesian network
bn <- motbf.fit(dag, data = trainingData, numIntervals = 4, POTENTIAL_TYPE = "MTE")

## -----------------------------------------------------------------------------
print(bn)

## -----------------------------------------------------------------------------
BiC.MoTBFBN(bn, data = testData)

## -----------------------------------------------------------------------------
## Obtain small training subset
set.seed(4)
dataTT <- TrainingandTestData(data, percentage_test = 0.99)

trainingData <- dataTT$Training
testData <- dataTT$Test
nrow(trainingData)

## -----------------------------------------------------------------------------
## Generate artificial prior dataset
means <- sapply(data, function(x){ifelse(is.numeric(x), mean(x),NA)})
set.seed(4)
priorData <- generateNormalPriorData(dag, data = trainingData, size = 5000, means = means)

## -----------------------------------------------------------------------------
## Learn univariate distribution using prior information
f <- learnMoTBFpriorInformation(priorData$aac, trainingData$aac, s = 5, POTENTIAL_TYPE = "MOP", returnAll = TRUE)

print(f)  

## -----------------------------------------------------------------------------
plot(f$posteriorFunction, xlim = f$domain, ylim = c(0,2.1), lwd = 3)
plot(f$dataFunction, xlim = f$domain, add = TRUE, col = 2, lwd = 3, lty = 2)
plot(f$priorFunction, xlim = f$domain, add = TRUE, col = 4, lwd = 3, lty = 3)
legend("topleft", legend = c("Posterior", "Data", "Prior"), col = c(1,2,4), lty = 1:3, lwd = 3, inset = c(0, -0.7), xpd = TRUE)

## -----------------------------------------------------------------------------
## Log-likelihood of the model that uses prior information
  sum(log(as.function(f$posteriorFunction)(testData$aac)))

## Log-likelihood of the model that does not use prior information
  sum(log(as.function(f$dataFunction)(testData$aac)))


## ----message=FALSE------------------------------------------------------------
## Fit Bayesian network using prior information
priorBN <- motbf.fit(dag, trainingData, numIntervals = 2,
                           POTENTIAL_TYPE = "MOP", s = 5, priorData = priorData)


## Fit BN without using prior information
BN <- motbf.fit(dag, trainingData, numIntervals = 2,
                           POTENTIAL_TYPE = "MOP")

# Compute log-likelihood 
logLikelihood.MoTBFBN(priorBN, data = testData)

logLikelihood.MoTBFBN(BN, data = testData)


## ----message=FALSE------------------------------------------------------------
## Learn a Bayesian network
dag <- LearningHC(data)

bn <- motbf.fit(dag, data = data, numIntervals = 4,
                         POTENTIAL_TYPE = "MOP")

## -----------------------------------------------------------------------------
# Specify the evidence set and target variable
obs <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 0.9, stringsAsFactors=FALSE)
node <- "alm2"

# Get the conditional distribution of 'node' and the generated sample
set.seed(4)
ap_post = get_approx_posterior(bn, target = node, 
                     evidence = obs, size = 100, maxParam = 8)

ap_post$fx

## ----message=F----------------------------------------------------------------
ex_post = variableElimination(bn, target = node, evidence = obs)


## ----messages=FALSE, warning=FALSE--------------------------------------------
# Plot the posterior distribution obtained with each solution
plot(ap_post$fx, col = "red") 
plot(ex_post, add = T, col = "blue")
legend("topleft", legend = c("Approximate", "Exact"), col = c("red", "blue"), lwd = 1,  inset = c(0, -0.5), xpd = TRUE)


## ----eval = requireNamespace("Rgraphviz", quietly = TRUE)---------------------
# Build a TAN model for classification
  bn_tan_cl = fit_tan("lip",data)
  bnlearn::graphviz.plot(getDAG(bn_tan_cl))


# Compute mutual information of each variable with mcg. 
# Firstly, discrete variables must be coerced to factor
  data[sapply(data, is.character)] = lapply(data[sapply(data, is.character)], as.factor)
  MI = mutual_information_tan(data,"mcg")
# Build a TAN model for regression
  bn_tan_reg = fit_tan("mcg",data, mutualInfoCond = MI)
  bnlearn::graphviz.plot(getDAG(bn_tan_reg))

## ----eval = requireNamespace("Rgraphviz", quietly = TRUE)---------------------
# Naive Bayes structure
nb = getStructure(data, "NB", "mcg")
bnlearn::graphviz.plot(nb)

# TAN wrapper (continuous variables are internally discretized)
tan_disc = getStructure(data, "TAN", "mcg")
bnlearn::graphviz.plot(tan_disc)

# HC wrapper
hc_bnlearn = getStructure(data, "HC")
bnlearn::graphviz.plot(hc_bnlearn)

