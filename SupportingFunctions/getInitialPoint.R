#' Initial Point
#' 
#' Finds an initial feasible point in the constraint set
#' 
#' @param index the indexing of the units in the experiment
#' @param Gamma the sensitivity parameter
#' 
#' @return a feasible rho and s pair
#' 
#' @export

getInitialPoint <- function(index, Gamma){
  # Find an initial feasible point in the constraint set
  
  I <- length(index)
  s <- rep(NA, I)
  rho <- list()
  for(i in 1:I){
    lix <- length(index[[i]])
    rho[[i]] <- rep(1, lix)*(1/lix)
    s[i] <- (1/lix + 1/(Gamma*lix))/2
  }
  list(rho=unlist(rho), s=s)
}
