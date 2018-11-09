Authors: Peter Cohen, Matt Olson, Colin Fogarty
Overview: The functions included in this project facilitate multivariate one-sided hypothesis tests for observational studies and conduct the associated sensitivity analyses.

Files:
  "chiBarSquaredTest.R" - a method to implement the chi-bar-squared test of Cohen, Olson, and Fogarty (2018)

  "chiBarSquaredTest_RealDataExample.R" - a commented example using "chiBarSquaredTest()" to perform a sensitivity analysis on real Data

  "chiBarSquaredTest_SyntheticDataExample.R" - a commented example using "chiBarSquaredTest()" to perform a sensitivity analysis on synthetic data (use can modify the synthetic data)

Subfolders:
  The "SupportingFunctions" folder contains:
    "calculateSigma.R" - a method for computing covariance matrices for a given set of unmeasured confounders

    "chibarsquare.R" - methods for various attributes of chi-bar-squared distribution (cdf, quantile function, weight function)

    "ChiBarSquareGradients.R" - methods for computing gradients of the chi-bar-squared cdf with respect to correlation matrix entries

    "computeTestStat.R" - a method to compute the chi-bar-squared test statistic of Cohen, Olson, and Fogarty (2018)

    "constrainProject.R" - a method to project infeasible solutions onto the feasible region (used in projected subgradient descent)

    "getInitialPoint.R" - a method to compute a feasible initial point for the projected subgradient descent

    "gradient.R" - a method to compute gradients for the projected subgradient descent

    "gradientDescent.R" - projected subgradient descent

    "innerSolve.R" - solves a quadratic program to compute optimal weights (lambdas) during projected subgradient descent

    "maxCorPair_ForPairedExperiments.R" and "minCorPair_ForPairedExperiments.R" - methods to compute the maximal (resp. minimal) correlation between two outcomes at a given Gamma in paired experiments

    "maxCorPair_FullMatching.R" and "minCorPair_FullMatching.R" - methods to compute the maximal (resp. minimal) correlation between two outcomes at a given Gamma in general fully matched experiments

    "multCompareFunctions.R" - methods from "http://www.mit.edu/~cfogarty/#software" based upon: Fogarty, C. and Small, D. (2016). Sensitivity Analysis for Multiple Comparisons in Matched Observational Studies through Quadratically Constrained Linear Programming. Journal of the American Statistical Association: Theory and Methods,111 (516), 1820-1830

    "permutationTest.R" - methods to perform a permutation test at Gamma = 1

  The "Synthetic Data Generation Functions" folder contains:
    "generateData3.R" - generates trivariate data at Gamma = 1 for a multitude of different experimental configurations
    "makeExperimentalSetupFromYData.R" - takes output from "generateData3()" and encodes the data in the form required for "chiBarSquaredTest()"
    "psiFunctions.R" - two psi functions for data processing
