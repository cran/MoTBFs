# MoTBFs
## Development version of the MoTBFs package

The aim of the **MoTBFs** package is to provide a free and accessible
implementation of algorithms for learning MoTBFs from data. The package
implements learning algorithms for univariate, conditional, and joint
MoTBF densities, explained in the [Mixtures of Truncated Basis Functions](vignettes/theoreticalBackgound.pdf) vignette. By extension, functionality is also provided for
learning MoTBF-based Bayesian networks by leveraging functionality from
the **bnlearn** package. Furthermore, the **MoTBFs** package also
provides methods for integrating prior domain knowledge in the learning
process, thus also supporting data sparse domains. An extensive example
included in the [Workflow of MoTBFs pacakge](vignettes/packageUsage.md) vignette illustrates the main functionality of the package.


New functionalities are included to the MoTBFs package, including the implementation of the variable elimination algorithm for inference, and the adaptation of the Chow-Liu algorithm to MoTBF densities to learn the structure of TAN models. 
