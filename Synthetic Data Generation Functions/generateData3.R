generateData3 <- function(rho, Tau_1, Tau_2, Tau_3, nostratum)
{
  errors = rmvnorm(n=nostratum, mean=c(0,0,0), sigma = rho + diag(x = 1 - rho, nrow = 3, ncol = 3))
  rootChiSquared = sqrt(rchisq(n = nostratum, df = 5) / 5)
  
  # Now, from Tau_1 and Tau_2 we can form the expected treated-minus-controol pair differences
  Ypairs_normal = data.frame(Y1 = errors[,1] + Tau_1,
                             Y2 = errors[,2] + Tau_2,
                             Y3 = errors[,3] + Tau_3)
  Ypairs_t = data.frame(Y1 = errors[,1]/rootChiSquared + sqrt(5 / 3)*Tau_1,
                        Y2 = errors[,2]/rootChiSquared + sqrt(5 / 3)*Tau_2,
                        Y3 = errors[,3]/rootChiSquared + sqrt(5 / 3)*Tau_3)
  
  s_k_normal = apply(abs(Ypairs_normal), MARGIN = 2, FUN = median) # s_k is set to the median currently
  
  #cat(s_k_normal)
  
  s_k_t = apply(abs(Ypairs_t), MARGIN = 2, FUN = median)  # s_k is set to the median currently
  
  normalData = list(YData = data.frame(Ypairs = Ypairs_normal), s_k = s_k_normal)
  tData = list(YData = data.frame(Ypairs = Ypairs_t), s_k = s_k_t)
  
  return(list(normalData = normalData, tData = tData)) 
}