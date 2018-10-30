################################################################################
#                     For Full-Matched experiments computes the maximal 
#                     correlation between two outcomes.
################################################################################

#INPUT:
#Q is a matrix with two outcomes
#index is an indexing set as elsewhere in this code
#Gamma is the sensitivity parameter

#Note: this performs the same task as maxCorPair_ForPairedExperiments.R but for general full-matchings
#       Since it is designed to run on general full-matchings, it runs slower than maxCorPair_ForPairedExperiments.R when #       run on pairs (there is less problem-specific structure to exploit).
#Note: the vector expyu is the vector with the (i,j)^th coordinate corresponding to exp(\gamma*u_{ij}).
maxCorFullMatching = function(Q, index, Gamma)
{

  epsilon = 1e-5
  K = 2 # We are finding pairwise correlations between exactly two of the outcome statistics
  nostratum = length(unique(index)) # Number of strata
  noIndiv = length(index) # Number of individuals
  
  CovPair_FullMatch = function(expyu, Q, index)
  {
    covariance = 0
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      qIthStrat = Q[ind, ]
      expyuIthStrat = expyu[ind] #the members of expyu that correspond to the individuals in the i^th stratum
      crossTerm = sum(qIthStrat[,1]*qIthStrat[,2]*expyuIthStrat) / sum(expyuIthStrat) #See page 1455 of Rosenbaum 2016
      productOfExpectations = sum(qIthStrat[,1]*expyuIthStrat) * sum(qIthStrat[,2]*expyuIthStrat) / (sum(expyuIthStrat)^2) #See page 1454 of Rosenbaum 2016
      covariance = covariance + (crossTerm - productOfExpectations)
    }
    return (covariance)
  }
  
  
  #computes the pairwise correlation
  CorPairMax_FullMatch = function(expyu, Q, index)
  {
    sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
    sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))
    corr.ret = CovPair_FullMatch(expyu, Q, index) / (sd1 * sd2)
    -corr.ret
  }
  

  #outputs the gradient of the covariance function with respect to the expyui terms
  GradCovPairMax_FullMatch = function(expyu, Q, index)
  {
    grad = rep(0, noIndiv)
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      qIthStrat = Q[ind, ]
      
      expyuIthStrat = expyu[ind]
      
      crossTerm = (qIthStrat[,1]*qIthStrat[,2]*sum(expyuIthStrat) - sum(qIthStrat[,1]*qIthStrat[,2]*expyuIthStrat)) / (sum(expyuIthStrat)^2)
      
      productTerm = ((qIthStrat[,1]*sum(qIthStrat[,2]*expyuIthStrat) + qIthStrat[,2]*sum(qIthStrat[,1]*expyuIthStrat)) * (sum(expyuIthStrat)^2) -
                       (sum(qIthStrat[,1]*expyuIthStrat)*sum(qIthStrat[,2]*expyuIthStrat)*2*sum(expyuIthStrat))) / (sum(expyuIthStrat)^4)
      
      grad[ind] = crossTerm - productTerm
    }
    return (grad)
  }
  
  #outputs the gradient of the pairwise correlation function with respect to the elements of expyu
  GradCorPairMax_FullMatch = function(expyu, Q, index)
  {
    sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
    sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))
    
    productDerivative = (sd2 * GradCovPairMax_FullMatch(expyu, cbind(Q[,1], Q[,1]), index) / (2 * sd1)) + (sd1 * GradCovPairMax_FullMatch(expyu, cbind(Q[,2], Q[,2]), index) / (2 * sd2))
    numerator = (GradCovPairMax_FullMatch(expyu, Q, index) * sd1 * sd2) - (CovPair_FullMatch(expyu, Q, index) * productDerivative)
    grad = numerator / ((sd1 * sd2)^2)
    -grad
  }
  

  expyu = runif(n = noIndiv, min = 1, max = Gamma) # random start point
  -optim(par = expyu, fn = CorPairMax_FullMatch, gr = GradCorPairMax_FullMatch, Q=Q, index=index, method = "L-BFGS-B", lower = 1 + epsilon, upper = Gamma - epsilon, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$value
}
