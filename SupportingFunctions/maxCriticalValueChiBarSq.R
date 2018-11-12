################################################################################
#       Computes the most conservative 1-alpha ChibarSq critical value at 
#       a given Gamma for the specified Q and index.
################################################################################
#index is given in the form of the matched set indexing

#' Maximal Chi-Bar-Sq. critical value
#' 
#' Computes the most conservative 1-alpha ChibarSq. critical value at a given Gamma for the specified Q and index. 
#' 
#' @param Q the matrix of the q_{ijk}.  It has K rows (the number of outcomes) and N columns (the number of individuals)
#' @param index an alternative form of indexing useful for computation.  The t^th element is the list of all individuals in the t^th matched set
#' @param Gamma the sensitivity parameter
#' @param alpha the significance level of the test
#' 
#' @return crit the most conservative feasible 1-alpha ChibarSq. critical value
#' 
#' @export


maxCritChiBarUB = function(Q, index, Gamma, alpha)
{
  K= ncol(Q)
  nostratum = length(unique(index))
  cmin = rep(0, K*(K-1)/2)
  cmax = cmin
  iter = 1
  
  pairMatching = (max(sort(table(index),decreasing=TRUE)) == 2) #indicator of whether the experiment is pair-matched
  

  if(pairMatching == FALSE)
  {
    for(j in 1:(K-1)) #The case of general (non-pairs) full matching.  For this we use the "FullMatching" methods
    {
      for(i in (j+1):K)
      {
        cmax[iter] = maxCorFullMatching(Q[,c(i,j)], index, Gamma) #the maximal correlation between the ith and jth outcomes
        cmin[iter] = minCorFullMatching(Q[,c(i,j)], index, Gamma) #the minimal correlation between the ith and jth outcomes
        iter = iter+1
      }
    }
  }else
  {
    for(j in 1:(K-1)) #The case of pairs matching.  For this we use the "pairs" methods (since they are faster)
    {
      for(i in (j+1):K)
      {
        cmax[iter] = maxCorPair(Q[,c(i,j)], index, Gamma) #the maximal correlation between the ith and jth outcomes
        cmin[iter] = minCorPair(Q[,c(i,j)], index, Gamma) #the minimal correlation between the ith and jth outcomes
        iter = iter+1
      }
    }
  }
  
  qtemp = qchisq(1-alpha, K)
  Cmin = diag(.5,K,K)
  Cmin[lower.tri(Cmin)] = cmin
  Cmin = Cmin + t(Cmin)
  
  Cmax = diag(.5,K,K)
  Cmax[lower.tri(Cmin)] = cmax
  Cmax = Cmax + t(Cmax)
  
  ret = 1
  CminInv = solve(Cmin)
  if(all(eigen(CminInv)$values >= 0))
  {
    qtemp = qchibarsq(1-alpha, CminInv)
  }
  
  #optimizer to find most conservative chibarsq critical value
  chibarcritopt = function(co, alpha, qtemp)
  {
    K = .5*sqrt(8*length(co) + 1) + .5
    Co = diag(.5,K,K)
    Co[lower.tri(Co)] = co
    Co = Co + t(Co)
    ret = 1
    CoInv = solve(Co)
    if(all(eigen(CoInv)$values >= 0))
    {
      if (!isTRUE(all.equal(CoInv, t(CoInv))) || any(diag(CoInv) < 0))
      {
        #Current CoInv is not a cov. matrix. (so we drive the optimization away from this parameter choice...somewhat like a Lagrangean approach to the constraint that the optimization only used covariance matrices)
        ret = 1
      }else{
        #Current CoInv is a cov. matrix.
        ret = pchibarsq(qtemp, CoInv)
      }
    }
    ret
  }
  
  crit = 0
  if(K == 2)
  {
    Sigworst = cmin + diag(1-cmin, K, K)
    crit = qchibarsq(1-alpha, solve(Sigworst))
  }else{
    cinit = .95*cmin+.05*cmax #starts very close to lower bound (experimentally we find that often the optimal solution is at the lower bound) this is just a vector of the lower triangular entries
    Cinit = .95*Cmin + .05*Cmax #actual matrix corresponding to cinit
    
    if((!(all(eigen(Cinit)$values >= 0))) && (!isTRUE(all.equal(Cinit, t(Cinit))) || any(diag(Cinit) < 0)))#if the Cinit picked is not a covariance matrix we use the cov matrix from Gamma = 1
    {
      nonMatchedSetIndexing = list() #converts from matched set indexing vector to the list of strata with corresponding inviduals in each list (this is called "index" elsewhere in the code, so the distinction here is made explicit)
      for(ind in 1:length(unique(index)))
      {
        nonMatchedSetIndexing[[ind]] = which(index == ind)
      }
      
      gamma1Point = getInitialPoint(index = nonMatchedSetIndexing, Gamma = 1)
      sigmaAtGamma1 = calculateSigma(rho = gamma1Point$rho, Q = t(Q), index = nonMatchedSetIndexing)
      cinit = as.vector(cov2cor(sigmaAtGamma1)[lower.tri(sigmaAtGamma1)])
    }
    
    
    cworst = cinit
    if(Gamma > 1) #optimizes over the feasible unmeasured confounders to find the most conservative critical value's associated inverse covariance matrix
    {
      cworst = optim(par = cinit, fn = chibarcritopt, alpha = alpha, qtemp = qtemp, method = "L-BFGS-B", lower = cmin, upper = cmax, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$par #no gradient usage (roughly the same performance as when gradients are used)
      
      # cworst = optim(par = cinit, fn = chibarcritopt, gr = chibarcritoptGradient, alpha = alpha, qtemp = qtemp, method = "L-BFGS-B", lower = cmin, upper = cmax, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$par #gradients used (included as a comment for completeness)
    }
    Cow = diag(.5,K,K)
    Cow[lower.tri(Cow)] = cworst
    Cow = Cow + t(Cow)
    ret = 1
    CowInv = solve(Cow)
    crit = qchibarsq(1-alpha, CowInv)
  }
  crit
}