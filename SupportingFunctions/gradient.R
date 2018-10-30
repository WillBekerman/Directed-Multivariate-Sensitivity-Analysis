gradient <- function(rho, lambda, Q, T, index){
  # Returns objective value and gradient of one-sided objective (given
  # optimal lambda!)
  # 
  # max_{lambda >= 0} (lambda'*(T-mu(rho))^2/lambda'*Sigma(rho)*lambda
  #
  # ARGS:
  #     rho:
  #     lambda:
  #     Q:
  #     T:
  #     index: list whose ith element contains vector of ith strata locations
  # OUTPUT:
  #       obj: objective value
  #       grad: [gradf; zeros(length(s),1]
  
  I <- length(index)
  N <- dim(Q)[2]
  
  # Use the following notation: h refers to numerator and g to the
  # denominator, so f = h/g, df = (gdh - hdg)/g^2
  alpha <- t(Q)%*%lambda
  g <- 0
  h <- (sum(lambda*T) - sum(alpha*rho))^2
  dh <- -2*alpha*(sum(lambda*T) - sum(alpha*rho))
  dg <- numeric(N)
  
  for(i in 1:I){
    ix <- index[[i]]
    rhoi <- rho[ix]
    alphai <- alpha[ix]
    g <- g + sum(rhoi*(alphai^2)) - sum(alphai*rhoi)^2
    dg[ix] <- alphai^2 - 2*alphai*sum(alphai*rhoi)
  }
  obj <- h/g
  grad <- (g*dh - h*dg)/(g^2)
  grad <- c(grad, rep(0, I))
  list(obj=obj, grad=grad)
  
}
