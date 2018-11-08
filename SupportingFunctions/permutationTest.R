################################################################################
#             Performs exact inference using permutation test (only for Gamma = 1)
################################################################################

#computes the test statistic for the Chibarsq test specifically at Gamma = 1
computeTestStatistic_Gamma1 = function(Q, TS, index, treatmentAllocation)
{
  EPS <- 1e-6

  initialPoint = getInitialPoint(index, Gamma = 1)
  
  rho = initialPoint$rho
  s = initialPoint$s
  

  #Set-up for optimization
  V <- calculateSigma(rho, Q, index)

  b <- (Q %*% treatmentAllocation) - (Q %*% rho) 
  K <- dim(V)[1]
  
  #solve inner optimization problem
  lambdaPos <- (solve.QP( Dmat = 2 * V,
                          dvec = 2 * b,
                          Amat = diag(1,K),
                          bvec = rep(0,K),
                          factorized = FALSE))$solution
  
  fpos <- sum(lambdaPos * b) / sqrt(as.vector(t(lambdaPos) %*% V %*% lambdaPos + EPS))
  
  return(fpos)
}

################################################################################
#                         Which block each indiv. belongs to
################################################################################
makeBlockIndices <- function(index)
{
  populationSize = max(index[[length(index)]])
  block.ind = unlist(apply(matrix(1:length(index)), MARGIN = 1, FUN = function(val)rep(val, length(index[[val]]))))
  
  return(block.ind)
}

################################################################################
#                         Performs a permutation test
################################################################################
permutationTest = function(Q, TS, index, direction = directionVector, alpha = alpha, Z=Z, subSampleSize = 500)
{
  populationSize = dim(Q)[2] #the number of individuals in the population (N in standard notation)

  ################################################################################
  #                           Calculate the optimal lambdas 
  #                           and observed test statistic
  ################################################################################
  observedTestStatistic = computeTestStatistic_Gamma1(Q, TS, index, Z)
  
  ################################################################################
  #                           Creates a collection of 
  #                           alternative treatment allocations 
  #                           (randomizes within strata)
  ################################################################################
  block.ind = makeBlockIndices(index) # Which block each indiv. belongs to
  B = length(unique(block.ind)) # The total number of blocks
  nt = tapply(Z, block.ind, sum) # The number of treated indiv. in each block
  nc = tapply(1-Z, block.ind, sum) # The number of control indiv. in each block
  
  randomizationMatrix = matrix(0, nrow = populationSize, ncol = subSampleSize)
  
  #creates a matrix of other randomized treatment allocations (respecting the initial stratification)
  for(s in 1:B)
  {
    tind = sample(which(block.ind==s), subSampleSize, replace = T)
    if(nt[s]==1)
    {
      randomizationMatrix[cbind(tind,1:subSampleSize)] = 1
    }else{
      randomizationMatrix[block.ind==s,] = 1
      randomizationMatrix[cbind(tind,1:subSampleSize)] = 0
    }
  }
  
  ################################################################################
  #                           Computes the test statistic w.r.t. 
  #                           the newly randomized treatment allocations
  ################################################################################
  randomizationOfTreatment = apply(randomizationMatrix, MARGIN = 2, FUN = function(val) computeTestStatistic_Gamma1(Q, TS, index, val))
  
  ################################################################################
  #                           Compute p-values from the permutation test
  ################################################################################
  alpha = .05 # The significance level
  pValue_Upper = (length(which(randomizationOfTreatment >= observedTestStatistic)) + 1) / (length(randomizationOfTreatment) + 1) #p-value for greater than test
  pValue_Lower = (length(which(randomizationOfTreatment <= observedTestStatistic)) + 1) / (length(randomizationOfTreatment) + 1) #p-value for less than test
  
  return(list(reject = (pValue_Upper <= alpha))) #currently just outputs the rejection for the greater-than test.
}