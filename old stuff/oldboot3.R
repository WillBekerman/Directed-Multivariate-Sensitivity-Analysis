# softmin <- function(x){ # gives min deviate
#   -log(sum(exp(-x)))
# }
# 
# softargmax <- function(x){ # gives lam
#   (exp(x)) / (sum(exp(x)))
# }
# 
# 
# EPS <- 1e-6
# V <- calculateSigma(rho, Q, index)
# b <- (TS - Q %*% rho)
# K <- dim(V)[1]
# 
# # Solves quadratic program over non-negative orthant
# lambdaVector <- solve.QP( Dmat = 2 * V,
#                             dvec = 2 * b,
#                             Amat = diag(1,K),
#                             bvec = rep(0,K))
# lambdaPos = lambdaVector$solution #solution over non-negative orthant
# 
# 
# 
# softargmax(softmin())
# 
# 
# 
# bprime = matrix(rep(seq(0,1,by=.0001),K),ncol=K)
# softmin(as.numeric(-2*t(b)%*%t(bprime)) + diag( bprime%*%V%*%t(bprime) ))
# 
# 
# 
# softargmax()






# Load necessary libraries
library(pracma)  # For numerical integration
library(nloptr)  # For numerical optimization
library(cubature)  # For multidimensional integration

# Define parameters
alpha <- 50  # Example value for alpha
t <- TS  # Example values for t (test statistics)
lambda_initial <- c(.01, .01, .01)  # Initial guess for lambda
N <- dim(Q)[2]
# Define the number of Monte Carlo samples
num_samples <- 10000

# Define mu(rho) and Sigma(rho) as functions of rho
mu <- function(rho) {
  # Ensure mu(rho) handles the dimensionality correctly
  return(Q %*% rho)
}

Sigma <- function(rho) {
  # Ensure Sigma(rho) handles the dimensionality correctly
  return(calculateSigma(rho, Q, index))
}

# Define f(rho)
f <- function(rho, lambda, t) {
  num <- sum(lambda * (t - mu(rho)))
  denom <- sqrt(sum(lambda * Sigma(rho) %*% lambda))
  return(num / denom)
}

# Define the derivative of exp(alpha * f(rho)) with respect to rho
derivative_exp_alpha_f <- function(rho, lambda, alpha, t) {
  f_val <- f(rho, lambda, t)
  exp_alpha_f <- exp(alpha * f_val)
  
  # Derivative of f(rho) with respect to lambda
  d_f_dlambda <- (t - mu(rho)) / sqrt(sum(lambda * Sigma(rho) %*% lambda)) - 
    (0.5 * f_val * Sigma(rho) %*% lambda / sum(lambda * Sigma(rho) %*% lambda))
  
  # Chain rule: derivative of exp(alpha * f(rho)) with respect to lambda
  d_exp_alpha_f_dlambda <- alpha * exp_alpha_f * d_f_dlambda
  
  return(d_exp_alpha_f_dlambda)
}

# # Integral of the derivative function
# integral_function <- function(lambda, alpha, t) {
#   integrand <- function(rho) {
#     derivative_exp_alpha_f(rho, lambda, alpha, t)
#   }
#   
#   # Multidimensional integration over the range of rho
#   result <- hcubature(integrand, lowerLimit = rep(0, N), upperLimit = rep(1, N))
#   
#   return(result$integral)
# }



# Monte Carlo approximation of the integral
monte_carlo_integral_function <- function(lambda, alpha, t, num_samples) {
  samples <- matrix(runif(num_samples * N), ncol = N)
  integral_sum <- 0
  for (i in 1:num_samples) {
    rho <- samples[i, ]
    integral_sum <- integral_sum + derivative_exp_alpha_f(rho, lambda, alpha, t)
  }
  return(integral_sum / num_samples)
}



# Objective function for optimization
objective_function <- function(lambda) {
  #abs(integral_function(lambda, alpha, t))
  abs(sum(monte_carlo_integral_function(lambda, alpha, t, num_samples)))
  
}

# Optimize lambda to minimize the objective function
result <- nloptr(
  x0 = lambda_initial,
  eval_f = objective_function,
  lb=rep(0,K),
  ub=rep(1,K),
  opts = list("algorithm" = "NLOPT_GN_ESCH", "xtol_rel" = 1e-8,"print_level"=3,maxeval=250)
)

# # Output the optimized lambda
# optimized_lambda <- result$solution
# 
# result <- nloptr(
#   x0 = optimized_lambda,
#   eval_f = objective_function,
#   lb=rep(0,K),
#   ub=rep(1,K),
#   opts = list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1e-8,"print_level"=3,maxeval=250)
# )
# # Output the optimized lambda
# optimized_lambda <- result$solution
print(optimized_lambda)

