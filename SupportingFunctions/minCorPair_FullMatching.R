################################################################################
#                     For Full-Matched experiments computes the minimal 
#                     correlation between two outcomes.
################################################################################

#INPUT:
#Q is a matrix with two outcomes
#index is an indexing set as elsewhere in this code
#Gamma is the sensitivity parameter

#Note: this performs the same task as minCorPair_ForPairedExperiments.R but for general full-matchings
#       Since it is designed to run on general full-matchings, it runs slower than minCorPair_ForPairedExperiments.R when #       run on pairs (there is less problem-specific structure to exploit).
#Note: the vector expyu is the vector with the (i,j)^th coordinate corresponding to exp(\gamma*u_{ij}).

#See maxCorPair_FullMatching.R for line-by-line comments

#' Computes maximal bivariate correlation
#' 
#' Computes the maximal correlation between two outcomes in a fully matched experiment
#' 
#' @param Q a subset of the original data matrix containing only two outcomes
#' @param index the indexing of the units in the experiment
#' @param Gamma the sensitivity parameter
#' 
#' @return the maximal correlation
#' @export

minCorFullMatching = function(Q, index, Gamma)
{
  epsilon = 1e-5
  K = 2
  nostratum = length(unique(index))
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
  CorPairMin_FullMatch = function(expyu, Q, index)
  {
    sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
    sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))
    corr.ret = CovPair_FullMatch(expyu, Q, index) / (sd1 * sd2)
    corr.ret
  }
  
  #outputs the gradient of the covariance function with respect to the expyui terms (updated)
  GradCovPairMin_FullMatch = function(expyu, Q, index)
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
  GradCorPairMin_FullMatch = function(expyu, Q, index)
  {
    sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
    sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))
    
    productDerivative = (sd2 * GradCovPairMin_FullMatch(expyu, cbind(Q[,1], Q[,1]), index) / (2 * sd1)) + (sd1 * GradCovPairMin_FullMatch(expyu, cbind(Q[,2], Q[,2]), index) / (2 * sd2))
    numerator = (GradCovPairMin_FullMatch(expyu, Q, index) * sd1 * sd2) - (CovPair_FullMatch(expyu, Q, index) * productDerivative)
    grad = numerator / ((sd1 * sd2)^2)
    grad
  }
  
  
  expyu = runif(n = noIndiv, min = 1, max = Gamma) # random start point
  optim(par = expyu, fn = CorPairMin_FullMatch, gr = GradCorPairMin_FullMatch, Q=Q, index=index, method = "L-BFGS-B", lower = 1 + epsilon, upper = Gamma - epsilon, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$value
}
