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