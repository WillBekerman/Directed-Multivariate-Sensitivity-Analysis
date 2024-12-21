
#----- fast evaluation of weights and quantiles of chibarsq dist. for known correlation -----

setwd("C:/Users/bekerman/Documents/Github/Directed-Multivariate-Sensitivity-Analysis")
source(file = "./chiBarSquaredTest.R")
library(tidyverse)
#sources all relevant files
allFiles = list.files(path = "./Synthetic Data Generation Functions")
for(file in allFiles)
{
fileName = paste("./Synthetic Data Generation Functions/", file, sep = "")
source(fileName)
}
allFiles = list.files(path = "./SupportingFunctions")
for(file in allFiles)
{
  fileName = paste("./SupportingFunctions/", file, sep = "")
  source(fileName)
}

alpha=1-.05
K=150
rho=0.2
nsim=1000000
componentvec=numeric(K+1)

for(sim in 1:nsim){

  #cat('\n\n\n\n\nSIMULATION NUMBER: ', sim, '\n\n\n\n\n\n')

  covmat = rho + diag(x = 1 - rho, nrow = K, ncol = K)
  y = fourPNO::rmvnorm(n=1, mu=rep(0,K), sigma = covmat)
  
  matprod=quadprog::solve.QP( Dmat = solve(covmat),
                    dvec = y%*%solve(covmat),
                    Amat = diag(1,K),
                    bvec = rep(0,K) )$solution
  
  numpos=sum(matprod>1E-6)
  componentvec[numpos+1] <- componentvec[numpos+1] + 1

}

wtsvec=componentvec/nsim
wts<-rev(wtsvec) #w_i(p,V) = w_{p-i}(p,V^{-1})
pchiroot = function(q, V, wts){pchibarsq(q, V, wts) - alpha}
ans = sqrt(uniroot(pchiroot, c(0, qchisq(alpha, K)), V=covmat, wts=wts)$root)

ans #sqrt(qchibarsq(alpha,solve(covmat))); 
wts #wchibarsq(solve(covmat))


