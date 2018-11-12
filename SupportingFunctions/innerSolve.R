# PURPOSE:
# solves max_{lambda} [lambda'*(T-mu(rho))]/sqrt(lambda'*Sigma(rho)*lambda)
#        s.t. lambda >= 0

# INPUT:
#       rho: the vector rho from Matt Olson's paper
#       Q: the matrix of q's from the test used
#       T: the observed test statistic
#       index: list of indices for who go treatment and who did not

# RETURNS:
#       lambda: optimal solution
#       fobj: optimal value


#' Solves inner optimization during subgradient descent
#' 
#' Uses dual method of Shapiro (2003) to solve the inner optimization for fixed rho
#' 
#' @param rho configuration of unmeasured confounders
#' @param Q the data matrix
#' @param TS the univariate test statistics
#' @param index the indexing of the units in the experiment
#' 
#' @return lambda: an optimal weighting 
#' @return  optval: the objective function
#' 
#' @export


innerSolve <- function(rho, Q, TS, index){
  EPS <- 1e-6
  V <- calculateSigma(rho, Q, index)
  b <- (TS - Q %*% rho)
  K <- dim(V)[1]
  
  # Solves quadratic program over non-negative orthant
  lambdaVector <- solve.QP( Dmat = 2 * V,
                               dvec = 2 * b,
                               Amat = diag(1,K),
                               bvec = rep(0,K))
  
  
  lambdaPos = lambdaVector$solution #solution over non-negative orthant
  UNconstrainedSoln = lambdaVector$unconstrained.solution #solution without any constraints

  
  fpos <- sum(lambdaPos * b) / sqrt(as.vector(t(lambdaPos) %*% V %*% lambdaPos + EPS)) #the objective value from the constrained solution
  
  if(fpos > 0)
  {
    lambda <- lambdaPos
    optval <- fpos 
  }else #in this case, we have found some confounder that puts our test satistic at 0 no matter what lambdas we pick. What this means is that we can stop with the optimization since we know that ther is some unmeasured confounder configuration that puts our test statitstic at the min possible value.
  {
    lambda <- rep(1, K)# Arbitrary choice
    optval <- 0
  }
  
  list(lambda=lambda, optval=optval)
}