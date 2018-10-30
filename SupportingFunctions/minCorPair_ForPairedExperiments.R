################################################################################
#                     For PAIRED experiments computes the minimal 
#                     correlation between two outcomes.
################################################################################

#INPUT:
#Q is a matrix with two outcomes
#index is an indexing set as elsewhere in this code, but must be pairs for this to output valid results
#Gamma is the sensitivity parameter

minCorPair = function(Q, index, Gamma)
{
  K = 2
  nostratum = length(unique(index))
  qnew = matrix(0, nostratum, K)
  for(j in 1:nostratum)
  {
    ind = which(index==j)
    qtemp = (Q[ind,])
    qnew[j,] = (qtemp[1,] - qtemp[2,])
  }
  
  CorPairMin = function(pprod, qnew, index)
  {
    qprod = qnew[,1]*qnew[,2]
    corr.ret = sum(pprod*qprod)/(sqrt(sum(pprod*qnew[,1]^2))*sqrt(sum(pprod*qnew[,2]^2)))
    corr.ret
  }
  
  
  GradCorPairMin = function(pprod, qnew, index)
  {
    
    qprod = qnew[,1]*qnew[,2]
    
    q1norm = sqrt(sum(pprod*qnew[,1]^2))
    q2norm = sqrt(sum(pprod*qnew[,2]^2))
    grad = (q1norm*q2norm*(qprod) - sum(qprod*pprod)*1/2*(q1norm/q2norm*qnew[,2]^2 + q2norm/q1norm*qnew[,1]^2))/(q1norm^2*q2norm^2)
    grad
  }
  
  
  
  pprod = rep((1/4+Gamma/(1+Gamma)^2)/2,nostratum)
  optim(pprod, CorPairMin, gr = GradCorPairMin, qnew=qnew, index=index, method = "L-BFGS-B", lower = Gamma/(1+Gamma)^2, upper = 1/4, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$value
}





