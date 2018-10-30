################################################################################
#       Computes the most conservative 1-alpha ChibarSq critical value at 
#       a given Gamma for the specified Q and index.
################################################################################
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
      ret = pchibarsq(qtemp, CoInv)
    }
    ret
  }
  
  crit = 0
  if(K == 2)
  {
    Sigworst = cmin + diag(1-cmin, K, K)
    crit = qchibarsq(1-alpha, solve(Sigworst))
  }else{
    cinit = .95*cmin+.05*cmax #starts very close to lower bound (experimentally we find that often the optimal solution is at the lower bound)
    cworst = cinit
    if(Gamma > 1) #optimizes over the feasible unmeasured confounders to find the most conservative critical value's associated inverse covariance matrix
    {
      cworst = optim(par = cmin, fn = chibarcritopt, alpha = alpha, qtemp = qtemp, method = "L-BFGS-B", lower = cmin, upper = cmax, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$par #no gradient usage (roughly the same performance as when gradients are used)
      
      # cworst = optim(par = cmin, fn = chibarcritopt, gr = chibarcritoptGradient, alpha = alpha, qtemp = qtemp, method = "L-BFGS-B", lower = cmin, upper = cmax, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$par #gradients used (included as a comment for completeness)
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