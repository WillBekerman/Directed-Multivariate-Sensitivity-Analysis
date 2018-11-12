#' Calculate covariance matrix 
#' 
#' Computes the covariance matrix for the given data at a fixed configuration of unmeasured confounders
#' 
#' @param Q the matrix of the q_{ijk}.  It has K rows (the number of outcomes) and N columns (the number of individuals)
#' @param rho the fixed set of unmeasured confounders
#' @param index an alternative form of indexing useful for computation.  The t^th element is the list of all individuals in the t^th matched set
#' 
#' @return Sigma: the covariance matrix
#' 
#' @export


calculateSigma <- function(rho, Q, index){
  # Calculates covariance matrix at a given rho
  
  I <- length(index) #number of strata
  K <- dim(Q)[1] #number of outcomes
  Sigma <- matrix(0, K, K) #initializes to zero matrix
  for(i in 1:I){
    ix <- index[[i]]
    rhoi <- rho[ix]
    Sigma <- Sigma + Q[,ix, drop=F] %*% diag(rhoi) %*% t(Q[,ix, drop=F])- 
      Q[,ix, drop=F] %*% outer(rhoi, rhoi) %*% t(Q[,ix, drop=F]) #see Cohen, Olson, and Fogarty for this equation
  }
  
  return(Sigma)
}