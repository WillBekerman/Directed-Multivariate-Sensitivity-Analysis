################################################################################
#                         Performs Optimization over lambda, rho, and s for
#                         Chi-Bar-Squared Test 
#                         (See appendix of Cohen, Olson, and Fogarty 2018)
################################################################################

#' Chi-Bar-Squared Test statistic
#' 
#' Performs Optimization over lambda, rho, and s for the Chi-Bar-Squared Test
#' 
#' @param Q the matrix of the q_{ijk}.  It has K rows (the number of outcomes) and N columns (the number of individuals)
#' @param TS the vector of the univariate test statistics
#' @param index an alternative form of indexing useful for computation.  The t^th element is the list of all individuals in the t^th matched set
#' @param Gamma the sensitivity parameter
#' @param Z the treatment indicator
#' @param alpha the significance level of the test (defaults to .05)
#' @param step the step-size in the subgradient descent (defaults to 100)
#' @param maxIter the maximum number of iterations of subgradient descent (defaults to 1000)
#' @param trueCrit the known, pre-computed chi-bar-squared critical value (defaults to NULL)
#' @param noCorBounds toggles optimistic speed-up using estimated worst-case rho to compute chi-bar-squared quantile if TRUE (defaults to FALSE)
#' @param useNormalQuantile toggles use of critical value from standard Normal distribution, instead of chi-bar-squared if TRUE (defaults to FALSE)

#' @return reject: indicator of rejection
#' @return lambdas: an optimal weighting of the outcomes
#' 
#' @export

computeTestStatistic = function(Q, TS, index, Gamma, Z, alpha, step, maxIter,
                                trueCrit, noCorBounds, useNormalQuantile)
{
  BETA <- 0.4 #Hyperparameter of optimization alg. (momentum term for gradient update)
  
  ################################################################################
  #                           Compute an initial feasible
  #                           soln. given index and Gamma
  ################################################################################
  outInit <- getInitialPoint(index, Gamma)
  
  ################################################################################
  #                           Perform Subgradient Descent
  ################################################################################
  outGrad <- gradientDescent(Q, TS, index, Gamma, rho0 = outInit$rho, s0 = outInit$s,
                             step = step, maxIter = maxIter, betam=BETA, alpha = alpha, Z=Z,
                             trueCrit = trueCrit, noCorBounds = noCorBounds, useNormalQuantile = useNormalQuantile)
  
  return(list(reject = outGrad$reject, contrast = outGrad$lambdas, invcovmat = outGrad$invcovmat))
}