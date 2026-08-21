Workflow of the MoTBFs package
================

The **MoTBFs** package is designed using *S3* objects. The functions
provided by the package implement the methods explained in the [Mixtures
of Truncated Basis Functions](theoreticalBackgound.html) vignette. The
package implements functions for learning univariate, multidimensional,
and conditional distributions, and provides support for parameter
learning in hybrid Bayesian networks. In addition, it includes functions
for incorporating prior knowledge when there is lack of data and for
carrying out probabilistic inference. Moreover, two classes are
incorporated in the package, `motbf` for defining univariate mixtures of
truncated basis functions and `jointmotbf` for specifying
multidimensional MoTBFs.

The functionality of the **MoTBFs** package is illustrated through an
analysis carried out on a real world dataset. More precisely, we use the
ecoli dataset \[1\], which is provided along with the package. The
dataset contains information about *Escherichia coli* and consists of
*n*=336 records, 8 input variables, and 1 output variable (the class).
It is a bacterium of the genus *Escherichia* that is commonly found in
the lower intestine of warm-blooded organisms. This dataset can be
downloaded from <http://archive.ics.uci.edu/ml/datasets/Ecoli>.

## How to install MoTBFs

The MoTBFs can be installed from CRAN, using the usual
`install.packages()` function.

``` r
## Install and load the MoTBFs package
 install.packages("MoTBFs")
 library("MoTBFs")

## Load ecoli dataset
 data("ecoli", package = "MoTBFs")
```

## The example dataset

The `ecoli` dataset is a data frame with 336 rows corresponding to
proteins and 9 columns corresponding to variables. The dataset contains
4 discrete variables, stored as characters, and 5 continuous variables.
The variables provide measurements of the cells used for predicting the
localization site of proteins. The first variable, `Sequence.Name`,
which is the accession number for the SWISS-PROT database, and the
output variable `class` will not be used in this running example, and we
will therefore remove them from the data frame. The discrete variables
`lip` and `chg` are binary attributes, where character numbers are used
as states; `"0.48"` and `"1"`, and `"0.5"` and `"1"`, respectively.

``` r
## Drop the first and last variables of the ecoli dataset 
  data <- ecoli[,-c(1,9)]
```

## Split the dataset into train and test

For validation purposes, the dataset is split into a training and a test
set using the `TrainingandTestData()` function.

``` r
## Split the dataset into train and test subsets
  set.seed(2)
  dataTT <- TrainingandTestData(data, percentage_test = 0.2)
  trainingData <- dataTT$Training
  testData <- dataTT$Test
```

The seed value determines the partitioning of the data into training and
test, and is therefore key to reproducing the experiments. From now on,
we will carry out all the analyses on the training data, leaving the
test dataset for estimating the predictive capabilities of the learned
models.

## Learn the directed acyclic graph

Our illustrative example basically consists of fitting MoTBF densities
to a previously learned Bayesian network structure over the variables in
the dataset. The structure can, for instance, be obtained, using the
function `hc()` from the **bnlearn** package. This function returns a
directed acyclic graph obtained from the dataset using a local search
method. For the sake of simplicity, we have included the function
`LearningHC()` in our package, which automatically converts into factors
those columns that are non-numeric, before calling the function `hc()`
in **bnlearn**. `LearningHC()` can also be used to discretize the
dataset before calling `hc()`, but we are not using this functionality
in the running example.

``` r
## Learn the structure of the Bayesian network using the training data
  dag <- LearningHC(trainingData)
  dag
#> 
#>   Bayesian network learned via Score-based methods
#> 
#>   model:
#>    [lip][alm1][mcg|lip:alm1][chg|lip][aac|alm1][gvh|mcg][alm2|gvh:lip:alm1] 
#>   nodes:                                 7 
#>   arcs:                                  8 
#>     undirected arcs:                     0 
#>     directed arcs:                       8 
#>   average markov blanket size:           3.14 
#>   average neighbourhood size:            2.29 
#>   average branching factor:              1.14 
#> 
#>   learning algorithm:                    Hill-Climbing 
#>   score:                                 BIC (cond. Gauss.) 
#>   penalization coefficient:              2.797356 
#>   tests used in the learning procedure:  102 
#>   optimized:                             TRUE
  
  # Plot the DAG using the graphviz.plot() from bnlearn package
  library(bnlearn)
  graphviz.plot(dag)
```

![](packageUsage_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Before describing how to learn the MoTBF distributions associated with
the network structure, we first present the basic functionality for
learning different types of MoTBF representations, i.e., univariate,
conditional, and joint MoTBF densities.

## Univariate MoTBFs densities

We illustrate the learning of a univariate MoTBF density by considering
the continuous variable `mcg`.

``` r
## Learn the density of variable mcg, using MTEs or MOP as basis functions
 f1 <- univMoTBF(trainingData[,1], POTENTIAL_TYPE = "MTE", nparam = 13)
 f2 <- univMoTBF(trainingData[,1], POTENTIAL_TYPE = "MOP", nparam = 11)
```

The `univMoTBF()` function is used for learning univariate densities.
The function is at the core of a collection of functions included in the
package to learn densities of class `motbf` from data. Least squares
optimization is used to minimize the mean squared error between the
empirical cumulative distribution and the estimated MoTBF.

The function takes two mandatory arguments, `data` and `POTENTIAL_TYPE`,
where the latter can either be `"MOP"` or `"MTE"` if polynomial or
exponential basis functions should be used, respectively. `univMoTBF()`
also accepts optional arguments: it is possible to specify the domain
over which the model will be fitted, `evalRange`, the exact number of
basis functions to be used, `nparam`, and the maximum number of
parameters in the function, `maxParam`, which selects the best fit using
the log-likelihood score. If `nparam` or `maxParam` are not given, then
the Bayesian information criterion (BIC) \[2\] is used for scoring and
function selection: it evaluates the two next functions and if the BIC
value does not improve then the function with the best BIC score so far
is returned.

The mathematical expression of the univariate density is shown via
`print()`, whereas other hidden elements related to the learning task
can be obtained using `summary()`.

``` r
 print(f2)
#> [1] 0.00528468532660471+30.7730282255632*x-1125.05786889453*x^2+16034.7596391062*x^3-118127.908988612*x^4+516328.233596937*x^5-1402763.10738281*x^6+2378343.85773264*x^7-2437579.93461867*x^8+1377857.84608826*x^9-329159.020877368*x^10
 summary(f2)
#> 
#>  MoTBFs FOR UNIVARIATE DISTRIBUTIONS 
#> 
#>  Model:
#>  0.00528468532660471+30.7730282255632*x-1125.05786889453*x^2+16034.7596391062*x^3-118127.908988612*x^4+516328.233596937*x^5-1402763.10738281*x^6+2378343.85773264*x^7-2437579.93461867*x^8+1377857.84608826*x^9-329159.020877368*x^10 
#> 
#>  Class: univmotbf mop motbf
#>  Subclass: mop 
#> 
#>  Coefficients:
#>  0.005284685 30.77303 -1125.058 16034.76 -118127.9 516328.2 -1402763 2378344 -2437580 1377858 -329159 
#> 
#>  Domain:
#>  (0, 0.89)
```

The object returned by `univMoTBF()` is a list of classes `univmotbf`, `motbf`, and either `mop` or `mte`, depending on the basis functions used. The object returned contais several
elements, including its mathematical expression and other hidden
elements related to the learning task. 

The learned densities can be plotted using the generic function
`plot()`. The next figure shows the the model fits, provided by
`univMoTBF()`, with blue dashed line for MOPs and red solid line for
MTEs overlaying the histogram of the training data of the `mcg`
variable.

``` r
## Plot the densities f1 and f2 over the histogram of variable mcg
  hist(trainingData[,1], prob = TRUE , main = "", xlab = "X", cex.lab=1.5, cex.axis=1.5)
  plot(f1, xlim = range(trainingData[,1]), col = "red", add = TRUE, lwd = 3)
  plot(f2, xlim = range(trainingData[,1]), col = "blue", add = TRUE, lwd = 3, lty = 2)
  legend("topleft", legend = c("MTE", "MOP"),col = c("red", "blue"), lty = 1:2, lwd = 3, cex = 1.75)
```

![](packageUsage_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

To evaluate the predictive ability of the models we use the generic
method `as.function()` developed for the `"motbf"` class to get the
log-likelihood as well as `BICMoTBF()` to obtain the BIC score.

``` r
## Compute log-likelihood and BIC score of the fitted densities
  sum(log(as.function(f1)(testData[,1])))
#> [1] 12.32783
  sum(log(as.function(f2)(testData[,1])))
#> [1] 10.51058
  
  BICMoTBF(f1,testData[,1])
#> [1] -17.10502
  BICMoTBF(f2,testData[,1])
#> [1] -14.71757
```

An alternative way to visually check the goodness of fit of the
estimated models is to simulate a data sample from the learned functions
and compare it with the training data. For doing this, we use the
inverse transform method, a technique for generating random samples from
a specific probability distribution based on evaluating the inverse of
the CDF on a uniform random number, yielding a value for the random
variable being sampled. This is done by function `rMoTBF()`. For the
sake of reproducibility, we fix the seed for the random numbers to be
used by the `rMoTBF()` function, which is set to 5 in this example. In
the next code snippet, the previous function fitted with a polynomial
basis, `f2`, will be used.

``` r
## Simulate data from the estimated density f2
  set.seed(5)
  X <- rMoTBF(size = 400, fx = f2)
  
## Test whether or not the simulated sample and the observed data come from the same distribution
  ks.test(trainingData[,1], X)
#> 
#>  Two-sample Kolmogorov-Smirnov test
#> 
#> data:  trainingData[, 1] and X
#> D = 0.065167, p-value = 0.5018
#> alternative hypothesis: two-sided
```

In this example the two-sample Kolmogorov-Smirnov test is used. The
*p*-value is notably above 0.05, so there is no evidence to reject the
null hypothesis that both samples are drawn from the same population.

<!-- We can compare the training data of variable `mcg` and the sample simulated from the distribution learned from the same data, using the `hist()` function. Moreover, we can compare the empirical cumulative distribution of the training and simulated data. -->

We can plot the histogram and the empirical cumulative distribution of
both the training data of variable `mcg` and the sample simulated from
the distribution learned from the same data, in order to compare both
distributions.

``` r
## Compare the histogram of both distributions
 hist(trainingData[,1], prob = TRUE, col = "#FFC107", density = 10, angle = -45)
 hist(X, prob = TRUE, col = "#0C7BDC", main = "", ylim = c(0,2.2), cex.lab=1.5, cex.axis=1.5, density = 10,add = T)
 legend('topleft', legend = c("Training data", "Simulated data"), fill = c("#FFC107", "#0C7BDC"), density = 20, angle = c(-45, 45))
```

![](packageUsage_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
## Compare the the CDF of both distributions
 plot(ecdf(trainingData[,1]), cex = 0, lwd = 3 , cex.lab = 1.5, cex.axis = 1.5, main = "")
 plot(integralMoTBF(f2), xlim = range(trainingData[,1]), col ="red", lwd = 3, add = TRUE)
 legend('topleft', legend = c("Training data", "Simulated data"), col = c("black", "red"), lwd = 3)
```

![](packageUsage_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

We can also manipulate the distributions with a collection of methods
for objects of class `"motbf"`. Here is an example of the use of three of them,
`coef()`, `integrate.motbf()`, and `derivMoTBF()`.

``` r
## Compute the derivative and integral of the fitted density

  # Coefficients of the fitted density
  coef(f2)
#>  [1]   5.284685e-03  3.077303e+01 -1.125058e+03  1.603476e+04 -1.181279e+05
#>  [6]   5.163282e+05 -1.402763e+06  2.378344e+06 -2.437580e+06  1.377858e+06
#> [11]  -3.291590e+05

  # Indefinite integral of the fitted density
  integrate.motbf(f2)
#> [1] 0.00528468532660471*x+15.3865141127816*x^2-375.01928963151*x^3+4008.68990977655*x^4-23625.5817977224*x^5+86054.7055994895*x^6-200394.729626116*x^7+297292.98221658*x^8-270842.21495763*x^9+137785.784608826*x^10-29923.547352488*x^11
  
  # Definite integral of the fitted density
  integralMoTBF(f2, lower = min(trainingData[,1]), upper = max(trainingData[,1]))
#> [1] 1
  
  # Derivarive of the fited density
  derivMoTBF(f2)
#> [1] 30.7730282255632-2250.11573778906*x+48104.2789173186*x^2-472511.635954448*x^3+2581641.16798468*x^4-8416578.64429686*x^5+16648407.0041285*x^6-19500639.4769494*x^7+12400720.6147943*x^8-3291590.20877368*x^9
```

## Joint MoTBFs densities

The learning process for multidimensional variables is similar to the
previous one. The function `jointmotbf.fit()` is used to solve the
quadratic optimization problem and returns the analytical expression of
the joint density. The returned object is of class `"jointmotbf"` and `"motbf"`.
The expression is the only visible element, while
the others can be retrieved using `attributes()`. In this example only
two variables are used, `mcg` and `alm1`, in order to be able to plot
the results.

``` r
## Learn joint distributions
  P = jointmotbf.fit(X = trainingData[,c("mcg", "alm1")], dimensions = c(5,5))
  
  attributes(P)
#> $names
#> [1] "Function"   "Domain"     "Iterations" "Time"      
#> 
#> $class
#> [1] "jointmotbf" "motbf"
```

The function `print()` can be used to obtain an expression of the
learned joint density, while `summary()` yields a more thorough excerpt
of the `"jointmotbf"` object.

``` r
 print(P)
#> [1] 1.00000000044876e-05-1.87890751245811e-13*alm1+1.34925514946196e-12*alm1^2-2.45341144148194e-12*alm1^3+1.31231637960609e-12*alm1^4-2.7069010790699*mcg+103.459764181587*mcg*alm1-461.084363497177*mcg*alm1^2+679.358798389865*mcg*alm1^3-319.267637943894*mcg*alm1^4-7.35251619834014*mcg^2+244.974957258816*mcg^2*alm1+38.6208995140072*mcg^2*alm1^2-1193.96736630721*mcg^2*alm1^3+920.72827509131*mcg^2*alm1^4+30.9216352421374*mcg^3-1102.42439690391*mcg^3*alm1+2410.7579267411*mcg^3*alm1^2-668.029611891877*mcg^3*alm1^3-677.36983290076*mcg^3*alm1^4-21.621357506081*mcg^4+782.648440776334*mcg^4*alm1-2103.42524400634*mcg^4*alm1^2+1294.26742623164*mcg^4*alm1^3+51.5825770400239*mcg^4*alm1^4

 summary(P)
#> 
#>  MoTBFs FOR MULTIVARIATE DISTRIBUTIONS 
#> 
#>  Model:
#>  1.00000000044876e-05-1.87890751245811e-13*alm1+1.34925514946196e-12*alm1^2-2.45341144148194e-12*alm1^3+1.31231637960609e-12*alm1^4-2.7069010790699*mcg+103.459764181587*mcg*alm1-461.084363497177*mcg*alm1^2+679.358798389865*mcg*alm1^3-319.267637943894*mcg*alm1^4-7.35251619834014*mcg^2+244.974957258816*mcg^2*alm1+38.6208995140072*mcg^2*alm1^2-1193.96736630721*mcg^2*alm1^3+920.72827509131*mcg^2*alm1^4+30.9216352421374*mcg^3-1102.42439690391*mcg^3*alm1+2410.7579267411*mcg^3*alm1^2-668.029611891877*mcg^3*alm1^3-677.36983290076*mcg^3*alm1^4-21.621357506081*mcg^4+782.648440776334*mcg^4*alm1-2103.42524400634*mcg^4*alm1^2+1294.26742623164*mcg^4*alm1^3+51.5825770400239*mcg^4*alm1^4
#> 
#>  Class: jointmotbf motbf
#> 
#>  Coefficients:
#>  1e-05 -1.878908e-13 1.349255e-12 -2.453411e-12 1.312316e-12 -2.706901 103.4598 -461.0844 679.3588 -319.2676 -7.352516 244.975 38.6209 -1193.967 920.7283 30.92164 -1102.424 2410.758 -668.0296 -677.3698 -21.62136 782.6484 -2103.425 1294.267 51.58258 
#> 
#>  Domain mcg:
#>  (0, 0.89)
#>  Domain alm1:
#>  (0.03, 1)
#> 
#>  Number of Iterations: 51 
#> 
#>  Processing Time: 0.005829096 secs
```

The processing time, `P$Time`, can vary depending on the CPU, 
but the learning outcome will always be the same for a specific data sample.

The generic function `plot()` can be used for objects of class `"jointmotbf"`.
This function accepts optional arguments such as `type`, where one can
choose between `"perspective"` and `"contour"`, `ranges`, used to
specify the plotting range, `orientation`, which indicates the
orientation of the perspective graph, and `filled` for getting a filled
contour plot.

``` r
## Plot the joint distribution of 2 variables
 par(mar=c(2,3,2,2))
 
 # Filled contour
 plot(P, data = trainingData[,c(1,6)])
 
 # Simple contour
 plot(P, data = trainingData[,c(1,6)], filled = FALSE, cex.lab = 2, cex.axis = 1.85, lwd = 1.5) 
  
 # Perspective
 plot(P, type = "perspective", data = trainingData[,c(1,6)], orientation=c(25,25), cex.lab = 2,xaxs = "i") 
```

<img src="packageUsage_files/figure-gfm/unnamed-chunk-16-1.png" width="50%" /><img src="packageUsage_files/figure-gfm/unnamed-chunk-16-2.png" width="50%" /><img src="packageUsage_files/figure-gfm/unnamed-chunk-16-3.png" width="50%" />

The `marginalJointMoTBF()` function computes the marginals of joint
densities. In this example we have two variables, so there are two
marginal densities.

``` r
## Compute the marginal distributions from the joint distribution (P)
 marginalJointMoTBF(P, var = 1)
#> [1] 0.00315094770414236+1.61192532573586*x+14.1692097860462*x^2-23.7487532213355*x^3+6.75404469435729*x^4
 marginalJointMoTBF(P, var = 2)
#> [1] -0.271699156295989+13.0463061825959*y-32.4520321021768*y^2+32.5558012281978*y^3-12.8781294161316*y^4
```

## Conditional MoTBFs densities

The next step in our analysis is learning conditional densities, which
is implemented by the function `conditionalMethod()`. Five of its
arguments are compulsory: `data`, the dataset; `nameParents`, a
character vector indicating the name of the parents; `nameChild`, a
character string containing the name of the child; `numIntervals`, the
maximum number of intervals for splitting the domain of the parent
variables; `POTENTIAL_TYPE`, the type of basis function. Other arguments
are optional, like `maxParam`, indicating the maximum number of
parameters for each function, and `s`, the expert’s relative confidence
in any prior knowledge, and `priorData` if prior knowledge is
incorporated in the analysis.

We will do the conditional analysis for only two variables in order to
be able to make a 2-dimensional plot of the obtained results. For
example, taking into account the relationship found by the dag, we
consider the child variable `gvh` with parent variable `mcg`.

``` r
## Learn conditional distributions
  P <- conditionalMethod(trainingData, nameParents = "mcg", nameChild = "gvh",
                         numIntervals = 5, POTENTIAL_TYPE ="MOP")
  printConditional(P)
#> Parent: mcg       Range: 0 < mcg < 0.44 
#> [1] 115.87045516994-1783.89439762351*x+10688.6519007417*x^2-32342.1750412397*x^3+54497.16895002*x^4-52033.2599582362*x^5+26393.6470625001*x^6-5536.00797133224*x^7
#> Parent: mcg       Range: 0.44 < mcg < 0.89 
#> [1] -37.300939128282+733.08773930547*x-5676.3802863816*x^2+22283.2874570132*x^3-47834.5916358261*x^4+56908.0376805823*x^5-35242.748222326*x^6+8867.55763320944*x^7
```

It can be noticed that the learning algorithm decides to split the
domain of the parent into two intervals even though we have set the
argument `numIntervals` to five. This is because the BIC score is not
improved any further by splitting the domain into more than two
intervals.

The resulting conditional density (a MOP in this case) can be plotted
using `plotConditional()`. The sample points can be overlaid by setting
the argument `points` to `TRUE`.

``` r
## Plot the conditional density of gvh given mcg
 plotConditional(P, data = trainingData, nameChild = "gvh", points = TRUE)
```

![](packageUsage_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## MoTBF distributions associated with the network structure

The last step is to learn the distributions tied to the Bayesian network
learned previously. For doing this task, the `MoTBFs_Learning()`
function of the **MoTBFs** package is used. The graph is a mandatory
argument, that can be of class `"bn"`, `"graphNEL"` or `"network"`.
Other mandatory arguments are the `data`, the maximum number of
intervals for splitting the domain of the parents (`numIntervals`), and
the type of basis function (`POTENTIAL_TYPE`). The function also accepts
additional arguments, but they are not listed here.

In the example, the DAG was obtained using the **bnlearn** package and
therefore it is an object of class `"bn"`. As an example, we will use a
maximum of 4 intervals and `"MTE"` potentials when learning the
densities (i.e. exponential basis functions).

``` r
## Learn the distributions of a Bayesian network
  bn <- MoTBFs_Learning(dag, data = trainingData, numIntervals = 4, POTENTIAL_TYPE = "MTE")
```

The results are reported using the `printBN()` function.

``` r
printBN(bn)
#> Potential(mcg)
#> Parent: alm1      Range: 0.03 < alm1 < 0.33 
#> Parent: lip       Range = "0.48" 
#> [1] 77.3522777593587-39.7786778378838*exp(2*x)-4.75251349022004*exp(-2*x)+7.48963890352964*exp(4*x)-110.446874605782*exp(-4*x)-0.485480169009668*exp(6*x)+71.340797051967*exp(-6*x)
#> Parent: alm1      Range: 0.33 < alm1 < 1 
#> Parent: lip       Range = "0.48" 
#> [1] -5.46201188285558+3.33472799427654*exp(2*x)+3.6655118677228*exp(-2*x)-0.423702303923044*exp(4*x)-1.08527005081554*exp(-4*x)
#> Parent: lip       Range = "1" 
#> [1] -11.4070647400056+6.00149378770896*exp(2*x)+4.44881770557646*exp(-2*x)-0.710921330310128*exp(4*x)+2.39443083395454*exp(-4*x)
#> 
#> Potential(gvh)
#> Parent: mcg       Range: 0 < mcg < 0.51 
#> [1] -97.6780018857865+34.1822579892428*exp(2*x)-282.216703494254*exp(-2*x)-2.72740482673279*exp(4*x)+2028.89862476474*exp(-4*x)-0.120964164514882*exp(6*x)-3609.50887689127*exp(-6*x)+0.0175533176036235*exp(8*x)+2069.44401943339*exp(-8*x)
#> Parent: mcg       Range: 0.51 < mcg < 0.89 
#> [1] 685.650778555292-360.289646374922*exp(2*x)+252.147763532602*exp(-2*x)+79.5432671120796*exp(4*x)-3074.91080772227*exp(-4*x)-8.35301855891214*exp(6*x)+4377.08778699041*exp(-6*x)+0.341014418584611*exp(8*x)-2004.85187262738*exp(-8*x)
#> 
#> Potential(lip)
#> 0.9630996 0.03690037 
#> 
#> 
#> Potential(chg)
#> Parent: lip       Range = "0.48" 
#> 1 0 
#> Parent: lip       Range = "1" 
#> 0.8181818 0.1818182 
#> 
#> 
#> Potential(aac)
#> Parent: alm1      Range: 0.03 < alm1 < 1 
#> [1] -3742.36654091738+2528.64291491454*exp(2*x)+1860.6477280111*exp(-2*x)-894.026826255328*exp(4*x)+2121.14505893388*exp(-4*x)+175.477131253675*exp(6*x)-3473.11261333144*exp(-6*x)-18.083692369602*exp(8*x)+1679.61261262166*exp(-8*x)+0.763490423609846*exp(10*x)-238.327761629406*exp(-10*x)
#> 
#> Potential(alm1)
#> [1] 158.612789592492-95.2652334081338*exp(2*x)-56.931839458767*exp(-2*x)+25.4631054109712*exp(4*x)-97.0824883095576*exp(-4*x)-3.162434543484*exp(6*x)+52.8052949749124*exp(-6*x)+0.14800462526047*exp(8*x)+17.628780248311*exp(-8*x)
#> 
#> Potential(alm2)
#> Parent: alm1      Range: 0.03 < alm1 < 0.33 
#> Parent: gvh       Range: 0.16 < gvh < 1 
#> Parent: lip       Range = "0.48" 
#> [1] 193.490320005721-112.656954852381*exp(2*x)-52.905855975255*exp(-2*x)+28.5291061114938*exp(4*x)-185.816773361928*exp(-4*x)-3.36983494270326*exp(6*x)+155.049264292973*exp(-6*x)+0.15173528062671*exp(8*x)-21.3863669285979*exp(-8*x)
#> Parent: alm1      Range: 0.33 < alm1 < 0.45 
#> Parent: gvh       Range: 0.16 < gvh < 1 
#> Parent: lip       Range = "0.48" 
#> [1] 395.060397879195-217.706142007792*exp(2*x)+16.0642632869484*exp(-2*x)+51.0248385142804*exp(4*x)-1009.41149259891*exp(-4*x)-5.62465460370298*exp(6*x)+1224.52541332394*exp(-6*x)+0.238712864505068*exp(8*x)-454.170336658462*exp(-8*x)
#> Parent: alm1      Range: 0.45 < alm1 < 0.71 
#> Parent: gvh       Range: 0.16 < gvh < 1 
#> Parent: lip       Range = "0.48" 
#> [1] -3.02605147368675+3.10977257300236*exp(2*x)+6.9062778763805*exp(-2*x)-0.801305572339328*exp(4*x)-12.57621579925*exp(-4*x)+0.0578364828307853*exp(6*x)+6.33068591306238*exp(-6*x)
#> Parent: lip       Range = "1" 
#> [1] 1.96482389541199-0.26607151573818*exp(2*x)-0.26607151573818*exp(-2*x)
#> Parent: alm1      Range: 0.71 < alm1 < 1 
#> Parent: gvh       Range: 0.16 < gvh < 1 
#> Parent: lip       Range = "0.48" 
#> [1] -1.26973277561207+0.635366387806036*exp(2*x)+0.635366387806036*exp(-2*x)
#> Parent: lip       Range = "1" 
#> [1] 1.01010101010101+0*exp(2*x)
```

Notice how nodes in the DAG with only discrete parents contain as many
functions as configurations of the parents, whereas nodes that have
continuous parents have at most 4 functions for each parent and,
finally, nodes that have mixed parents contain as many functions as
configurations of the discrete parents times the number of regions into
which the domain of the continuous parents is split.

The BIC criterion is used to decide the number of splitting points of
the domain of the continuous parent nodes and to choose the number of
basis functions used. The function `BiC.MoTBFBN()` can be used to
compute the log-likelihood and the BIC score of a dataset given the
Bayesian network.

``` r
BiC.MoTBFBN(bn, data = testData)
#> $LogLikelihood
#> [1] 173.5496
#> 
#> $BIC
#> [1] -51.40147
```

### Learn MoTBF distributions in a full hybrid network using prior knowledge

We will now exemplify the use of prior knowledge in the learning
process. In order to illustrate the approach, we first select a small
subset of the `Ecoli` dataset using `TrainingandTestData()`. In the next
example the percentage of the test data is 99%, which means the training
data is only 1% of the full dataset.

``` r
## Obtain small training subset
  set.seed(4)
  dataTT <- TrainingandTestData(data, percentage_test = 0.99)
  
  trainingData <- dataTT$Training
  testData <- dataTT$Test
  nrow(trainingData)
#> [1] 13
```

There are 13 entries in the training dataset. We are going to fit MoTBFs
with and without prior information.

Learning univariate and conditional distributions and Bayesian networks
can be done using the functions `learnMoTBFpriorInformation()` and
`MoTBFs\_Learning()`. The arguments for these functions are the same as
previously explained and, in addition, it is necessary to specify the
expert confidence in the prior knowledge, `s`, and the prior dataset
`priorData`. On the one hand, to generate an artificial prior dataset
the `generateNormalPriorData()` function can be used.

``` r
## Generate artificial prior dataset
  means <- sapply(data, mean)
  set.seed(4)
  priorData <- generateNormalPriorData(dag, data = trainingData, size = 5000, means = means)
```

On the other hand, argument *s* takes values on the interval \[0,*N*\],
where *N* is the sample size, and is used to synchronize the support of
the prior knowledge and the sample. We refer the reader to \[3\] for the
details. In this example we will use the `aac` variable from the data
set, have `s = 5` as confidence level, and set `"MOP"` as potential
type.

``` r
## Learn univariate distribution using prior information
  f <- learnMoTBFpriorInformation(priorData$aac, trainingData$aac, s = 5, POTENTIAL_TYPE = "MOP")
  
  print(f)  
#> $coeffs
#> [1] 0.5509206 0.4490794
#> 
#> $posteriorFunction
#> [1] 0.2911610116964+2.76267070095529*x+4.37213811054131*x^2-30.4218251545051*x^3-0.838982034246077*x^4+965.84448245364*x^5-4183.3860786611*x^6+7735.84337198152*x^7-7321.9450341575*x^8+3498.94490221637*x^9-671.510406730835*x^10
#> 
#> $priorFunction
#> [1] 0.0673759470400333+1.76398164857178*x+11.8532729256404*x^2-55.2199852318976*x^3-1.52287298035552*x^4+1753.14655799018*x^5-7593.44701738851*x^6+14041.6676050034*x^7-13290.3826315979*x^8+6351.08790634143*x^9-1218.88790545662*x^10
#> 
#> $dataFunction
#> [1] 0.565695496832341+3.9878399732122*x-4.80554985469248*x^2
#> 
#> $domain
#> [1] -0.1232877  0.9531283
```

`learnMoTBFpriorInformation()` returns the fitted model using the
training data only (`$dataFunction`), the model fit using the prior data
only (`$priorFunction`) and the model fit that combines the training
data and the prior data (`$posteriorFunction`). The three univariate
densities can be plotted using the generic method `plot()`.

``` r
  plot(f$posteriorFunction, xlim = f$domain, ylim = c(0,2.1), lwd = 3)
  plot(f$dataFunction, xlim = f$domain, add = TRUE, col = 2, lwd = 3, lty = 2)
  plot(f$priorFunction, xlim = f$domain, add = TRUE, col = 4, lwd = 3, lty = 3)
  legend("topleft", legend = c("Posterior", "Data", "Prior"), col = c(1,2,4), lty = 1:3, lwd = 3)
```

![](packageUsage_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
## Log-likelihood of the model that uses prior information
  sum(log(as.function(f$posteriorFunction)(testData$aac)))
#> [1] 134.566
## Log-likelihood of the model that does not use prior information
  sum(log(as.function(f$dataFunction)(testData$aac)))
#> [1] 78.96405
```

The best model, taking into account the log-likelihood, is the MoTBF
which uses the prior data, `f$posteriorFunction`.

The last step is to incorporate the prior knowledge in the full Bayesian
network. For this analysis we are not going to print out the results
(which could be done using function `printBN()`), because the structure
is similar to the previous Bayesian network representations. As an
example, we will use `numIntervals = 2`, `POTENTIAL_TYPE = "MOP"`, and
`s = 5`.

``` r
## Fit Bayesian network using prior information
priorBN <- MoTBFs_Learning(dag, trainingData, numIntervals = 2,
                           POTENTIAL_TYPE = "MOP", s = 5, priorData = priorData)


## Fit BN without using prior information
BN <- MoTBFs_Learning(dag, trainingData, numIntervals = 2,
                           POTENTIAL_TYPE = "MOP")

# Compute log-likelihood 
logLikelihood.MoTBFBN(priorBN, data = testData)
#> [1] 124.384
logLikelihood.MoTBFBN(BN, data = testData)
#> [1] 14.64589
```

Looking at the log-likelihood corresponding to the network with and
without prior data, we can see that, in this example, incorporating
prior knowledge is better when data is scarce.

## Inference

After a Bayesian network has been constructed, the **MoTBFs** package
can be used to obtain the conditional density of any variable in the
network given that some other variables have been observed. The
conditional distribution is obtained by forward sampling. As an example,
consider a network estimated from the `ecoli` dataset:

``` r
## Learn a Bayesian network
  dag <- LearningHC(data)
  
  bn <- MoTBFs_Learning(dag, data = data, numIntervals = 4,
                           POTENTIAL_TYPE = "MTE")
  
```

The observed values are specified using a data frame. In the example, we
are assuming that we want to compute the conditional density of `alm2`
given that `lip="0.48"`, `alm1 = 0.55` and `gvh = 1`. This is achieved
by using the function `forward_sampling()` were we have chosen a sample
size equal to 10 specified by parameter `size = 10`.

``` r
## Specify the evidence set and target variable
  obs <- data.frame(lip = "0.48", alm1 = 0.55, gvh = 1, stringsAsFactors=FALSE)
  node <- "alm2"

## Get the conditional distribution of 'node' and the generated sample
  set.seed(5)
  forward_sampling(bn, dag, target = node, evi = obs, size = 10, maxParam = 15)
#> Processing Time: 0.0895860195159912secs
#> $fx
#> [1] -4.37388075852153+2.45524553272302*exp(2*x)+1.50544468385721*exp(-2*x)-0.23926467845802*exp(4*x)+2.00453598788814*exp(-4*x)
#> 
#> $sample
#>          mcg gvh  lip chg       aac alm1       alm2
#> 1  0.7450156   1 0.48 0.5 0.3564026 0.55 0.53709571
#> 2  0.2493266   1 0.48 0.5 0.4873396 0.55 0.55395502
#> 3  0.4075408   1 0.48 0.5 0.5052861 0.55 0.74832025
#> 4  0.3205169   1 0.48 0.5 0.4844040 0.55 0.82546475
#> 5  0.4448844   1 0.48 0.5 0.4248584 0.55 0.07672959
#> 6  0.5717975   1 0.48 0.5 0.4611889 0.55 0.35412707
#> 7  0.7548104   1 0.48 0.5 0.5695499 0.55 0.68322911
#> 8  0.5475842   1 0.48 0.5 0.7797320 0.55 0.45343460
#> 9  0.3048408   1 0.48 0.5 0.5251502 0.55 0.73064241
#> 10 0.7281756   1 0.48 0.5 0.3819332 0.55 0.70396470
```

The output consists of the posterior density and the sample from which
the density parameters were estimated.

## References

<div id="refs" class="references csl-bib-body" line-spacing="2">

<div id="ref-Lichman:2013" class="csl-entry">

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Lichman, M. (2013). UCI machine learning
repository. University of California, Irvine, School of Information;
Computer Sciences. Retrieved from <http://archive.ics.uci.edu/ml></span>

</div>

<div id="ref-Sch78" class="csl-entry">

<span class="csl-left-margin">2. </span><span
class="csl-right-inline">Schwarz, G. (1978). Estimating the dimension of
a model. *Annals of Statistics*, *6*, 461–464.</span>

</div>

<div id="ref-Per15" class="csl-entry">

<span class="csl-left-margin">3. </span><span
class="csl-right-inline">Pérez-Bernabé, I., Fernández, A., Rumí, R., &
Salmerón, A. (2016). Parameter learning in hybrid Bayesian networks
using prior knowledge. *Data Mining and Knowledge Discovery*, *30*,
576–604.</span>

</div>

</div>
