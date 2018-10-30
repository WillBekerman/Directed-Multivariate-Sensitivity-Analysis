################################################################################
#       Gradients for Chi-Bar-Squared CDF as a function of correlations
################################################################################
# For the relevent paper, see "A Multivariate Analogue of the One-Sided Test" by Akio Kudo (1963)
# Compute Gradients of Chi-bar-squared Weights with respect to correlation matrix entries
# INPUTS:
  # k: the number of 'degrees of freedom' of the Chi-bar-squared distribution
  # Correlation: the correlation structure (V in the notation of Sen and Silvapulle)

chibarcritoptGradient = function(co, alpha, qtemp)
{
  K = .5*sqrt(8*length(co) + 1) + .5
  Co = diag(.5,K,K)
  Co[lower.tri(Co)] = co
  Co = Co + t(Co)
  ret = 1
  CoInv = solve(Co)
  if(all(eigen(CoInv)$values >= 0))
  {
    ret = pchibarsqGradient(qtemp, CoInv)
    ret = as.vector(ret[lower.tri(ret)])
  }
  ret
}

################################################################################
#Gradients for the Chi-Bar-Squared CDF as a function of the correlations Helper Methods
################################################################################
# INPUTS:
  # 'q' is the argument of the CDF
  # 'V' is the covariance matrix
  # 'lower.tail' toggles the lower tail (CDF) or upper tail (1 - CDF)

# RETURNS:
  # 'grad' the gradient vector with respect to the correlations
pchibarsqGradient = function(q, V, lower.tail = TRUE)
{
  k = nrow(V)
  
  grad = pchisq(q, 0, lower.tail=FALSE)*iThGradient(0, k, V)
  
  for(i in 1L:k){
    grad = grad + pchisq(q, i, lower.tail=FALSE)*iThGradient(i, k, V)
  }
  
  grad = if(isTRUE(lower.tail)) (-1)*grad else grad
  
  return(grad)
}

#Compute gradient of i^th weight w_{i}(k, V)
iThGradient = function(i, k, V)
{
  stopifnot(i >= 0 && i <= k) #safety check
  
  if (i > 0 && i < k)
  {
    V = cov2cor(V) #just to be safe in case other methods require correlation matrices
    subsets = combn(1:k, i) #lists all length i subsets of {1, ..., k}. Each one is its own column
    
    grad = rowSums(apply(subsets, MARGIN = 2, FUN = function(x){iThGradient_GivenM(i, k, V, x)}))
    
    #depending on formatting issues it might be better to output this as a matrix
    grad = matrix(grad, nrow = k, ncol = k)
    return(grad)
  }else if (i == 0)
  {
    N = 1:k
    M = c() #empty set
    LambdaN = V[N, N]
    LambdaN_inv = solve(LambdaN)
    grad = gradientForLambda_N_Inverse(i, k, V, M)
    return(grad)
  }else #i = k case
  {
    M = 1:k
    N = c() #empty set
    grad = gradientForLambda_MN(i, k, V, M)
    return(grad)
  }
}


#computes the gradients of the term with a fixed given 'M' in the sum over all subsets {1, ..., k}
iThGradient_GivenM = function(i, k, V, M)
{
  N = setdiff(1:k, M) #N = {1, ..., k} \setminus M
  
  LambdaN = V[N, N]
  LambdaN_inv = solve(LambdaN)
  LambdaMN = V[M, M] - V[M, N] %*% solve(LambdaN) %*% V[N, M] # Check this
  
  grad = gradientForLambda_N_Inverse(i, k, V, M) * pmvnorm(lower = rep(0, i), sigma = LambdaMN) +
    pmvnorm(lower = rep(0, k - i), sigma = LambdaN_inv)*gradientForLambda_MN(i, k, V, M)
  return(grad)
}

LambdaNTerm = function(i, k, V, M)
{
  N = setdiff(1:k, M) #N = {1, ..., k} \setminus M
  
  LambdaN = V[N, N]
  
  LambdaN_inv = solve(LambdaN) #pre-computed to save time later
  LambdaN_inv_corr = cov2cor(LambdaN_inv) #pre-computed to save time later
  
  return (pmvnorm(lower = rep(0, k - i), sigma = LambdaN_inv_corr))
}

gradientForLambda_N_Inverse = function(i, k, V, M)
{
  N = setdiff(1:k, M) #N = {1, ..., k} \setminus M
  
  LambdaN = V[N, N]
  
  LambdaN_inv = solve(LambdaN) #pre-computed to save time later
  LambdaN_inv_corr = cov2cor(LambdaN_inv) #pre-computed to save time later
  
  # LambdaMN = V[M, M] - V[M, N] %*% LambdaN_inv %*% V[N, M] # Not needed in this method
  
  grad = 0 * V # Just a convenient way to get a zero matrix of the right dimensions
  
  JacobianOfPDF = pmvnorm(lower = rep(0, k - i), sigma = LambdaN_inv_corr) * (-1) * as.vector(LambdaN)
  
  for (s in 1:k)
  {
    for (t in 1:k)
    {
      Mbar = diag(x = 0, nrow = length(N)) # Get a matrix of all zeros of the correct size
      if (s %in% N && t %in% N)
      {
        i = which(N == s)
        j = which(N == t)
        Mbar[i, j] = 1
        Mbar[j, i] = 1
      }
      JacobianOfInverse = as.vector(-LambdaN_inv %*% Mbar %*% LambdaN_inv)
      grad[s, t] = JacobianOfPDF %*% JacobianOfInverse #See eqn (5) in the notes
    }
  }
  
  return(grad)
}

gradientForLambda_MN = function(i, k, V, M)
{
  if (i < k)
  {
    N = setdiff(1:k, M) #N = {1, ..., k} \setminus M
    
    LambdaN = V[N, N]
    LambdaN_inv = solve(LambdaN)
    LambdaMN = V[M, M] - V[M, N] %*% solve(LambdaN) %*% V[N, M] # Check this
  }else
  {
    LambdaMN = V[M, M] #this case has no N set to condition on
  }
  
  
  grad = 0 * V # Just a convenient way to get a zero matrix of the right dimensions
  
  JacobianOfPDF = pmvnorm(lower = rep(0, i), sigma = LambdaMN) * (-1) * as.vector(solve(LambdaMN))
  
  for (s in 1:k)
  {
    for (t in 1:k)
    {
      if (s %in% M && t %in% M) #s,t in M case, notice that this is the only case that matters when M = 1:k
      {
        Mbar = diag(x = 0, nrow = length(M)) # Get a matrix of all zeros of the correct size
        i = which(M == s)
        j = which(M == t)
        Mbar[i, j] = 1
        Mbar[j, i] = 1
        JacobianOfSecondTerm = as.vector(Mbar)
        grad[s, t] = JacobianOfPDF %*% JacobianOfSecondTerm #See eqn (5) in the notes
      }else if (s %in% N && t %in% N) #s,t in N case
      {
        Mbar = diag(x = 0, nrow = length(N)) # Get a matrix of all zeros of the correct size
        i = which(N == s)
        j = which(N == t)
        Mbar[i, j] = 1
        Mbar[j, i] = 1
        JacobianOfSecondTerm = as.vector(V[M, N] %*% LambdaN_inv %*% Mbar %*% LambdaN_inv %*% V[N, M])
        grad[s, t] = JacobianOfPDF %*% JacobianOfSecondTerm #See eqn (5) in the notes
      }else #(s in M and t in N) or (s in N and t in M) cases
      {
        Mbar = matrix(0, nrow = length(M), ncol = length(N))
        if (s %in% M && t %in% N)
        {
          i = which(M == s)
          j = which(N == t)
        }
        if (s %in% N && t %in% M)
        {
          i = which(M == t)
          j = which(N == s)
        }
        Mbar[i, j] = 1
        JacobianOfSecondTerm = -as.vector(Mbar %*% LambdaN_inv %*% V[N, M] + V[M, N] %*% LambdaN_inv %*% t(Mbar))
        grad[s, t] = JacobianOfPDF %*% JacobianOfSecondTerm #See eqn (5) in the notes
      }
    }
  }
  
  return(grad)
}
