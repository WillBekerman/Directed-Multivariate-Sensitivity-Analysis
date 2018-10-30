################################################################################
#                         Performs Optimization over lambda, rho, and s for
#                         Chi-Bar-Squared Test 
#                         (See appendix of Cohen, Olson, and Fogarty 2018)
################################################################################
computeTestStatistic = function(Q, TS, index, Gamma, direction, Z,
                                             alpha,  step, maxIter)
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
                            step = step, maxIter = maxIter, betam=BETA, alpha = alpha, Z=Z)
  
  return(list(reject = outGrad$reject, lambdas = outGrad$lambdas))
}