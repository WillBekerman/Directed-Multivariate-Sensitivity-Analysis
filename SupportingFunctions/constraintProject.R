# solves min_{x,s} ||(x,s) - g||^2
#       s.t. 0 <= x <= 1, sum(x) = 1, s <= x_i <= Gamma*s
#
# In practice, g will be a gradient step of the form 
# (x, s) - t*(dx_J, ds_j)
#
# ARGS:
#     g: gradient step of the form (x, s) - t*(dx_J, ds_j)
#
# OUTPUT:
#       out: list with the projection, elements s and rho
#
constraintProject <- function(Gamma, g){
  n <- length(g) - 1
  Dmat <- diag(rep(1,n+1))
  dvec <- g
  bvec <- c(1, rep(0, n), rep(-1, n), rep(0, n), rep(0,n))
  Amat <- cbind(c(rep(1,n), 0),
                rbind(diag(1,n), 0),
                rbind(-diag(1,n), 0),
                rbind(diag(1,n), -1),
                rbind(-diag(1,n), Gamma))
  out <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  list(s = out$solution[n+1], x= out$solution[1:n])
}

