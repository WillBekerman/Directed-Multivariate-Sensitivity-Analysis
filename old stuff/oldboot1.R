set.seed(1)
source(file = "./chiBarSquaredTest.R")
library(tidyverse)

#sources all files for generating synthetic data
allFiles = list.files(path = "./Synthetic Data Generation Functions")
for(file in allFiles)
{
  fileName = paste("./Synthetic Data Generation Functions/", file, sep = "")
  source(fileName)
}

################################################################################
#                             Tweakable Parameters
################################################################################
#Total number of Gamma's to try with the double-and-halve approach
numGamma = 10

#Epsilon (if the next Gamma differs from the current Gamma by < epsilon then quit)
epsilon = 1E-3

#The significance level of the test
alpha = .05

#Toggles printout amount
verbose = T

#Toggles diagnostics
showDiagnostics = T

#number of strata
nostratum = 500#5000
# 
#tau_1 (impact on first outcome)
Tau_1 = -.75

#tau_2 (impact on second outcome)
Tau_2 = .5

#tau_3 (impact on third outcome)
Tau_3 = .75

#Taus <- c(-0.5, -0.5, 0.25, 0.25, 0.5, 0.5)

#correlation
correlation = 0

#directions
directions = c("Less", "Greater", "Greater")
#directions=c(rep("Less",2),rep("Greater",4))

#split proportion
planning_sample_prop=0.2

#divide data into planning and analysis samples
planning_sample_size_sets <- floor(planning_sample_prop*nostratum)
ix_planning <- sort(sample(x=1:nostratum,size=planning_sample_size_sets))
ix_analysis <- (1:nostratum)[-ix_planning]
n1 <- planning_sample_size_sets
n2 <- nostratum-planning_sample_size_sets








syntheticData = generateData3(rho = correlation, Tau_1 = Tau_1, Tau_2 = Tau_2, Tau_3 = Tau_3, nostratum = nostratum)





nboot=1000
lamlist=vector('list',nboot)
dessensvec=numeric(nboot)
for (b in 1:nboot){
  Data_asymp = syntheticData$normalData
  Data_asymp$YData <- Data_asymp$YData[ix_planning,]
  
  
  
  bootix=sample(1:nrow(Data_asymp$YData),size=nrow(Data_asymp$YData),replace=T)
  bootdat=Data_asymp$YData[bootix,]
  Data_asymp$YData=bootdat
  
  Data_asymp$s_k = apply(abs(Data_asymp$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
  Data_asymp$omega_k = Data_asymp$s_k  # omega_k is set to the median currently
  psi_asymp = psi_huber #this can be changed to use an inner-trimmed psi function
  experimentalSetup_asymp = makeExperimentalSetupFromYData_asymp(Data_asymp, psi_asymp)
  psi_YOverOmega_asymp = experimentalSetup_asymp$psi_YOverOmega
  psi_YOverOmega_asymp[,directions == 'Less']=-psi_YOverOmega_asymp[,directions == 'Less']
  K_asymp=ncol(psi_YOverOmega_asymp)
  eta <- function(lam,psi_res){
    mean(sqrt( (lam%*%t(psi_res))^2 ))
  }
  theta <- function(lam,psi_res){
    mean( lam%*%t(psi_res) )
  }
  get_design_sens <- function(lam,psi_res){
    (eta(lam,psi_res)+theta(lam,psi_res)) / (eta(lam,psi_res)-theta(lam,psi_res))
  }
  lam_optdesignsens <- try({
    optim(par = rep(0.01,K_asymp), fn = get_design_sens, psi_res = psi_YOverOmega_asymp, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1e-3))$par
  })
  if (class(lam_optdesignsens) == "try-error") {
    lam_optdesignsens <- try({
      optim(par = rep(0.01,K_asymp), fn = get_design_sens, psi_res = psi_YOverOmega_asymp, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1))$par})
  }
  if (class(lam_optdesignsens) == "try-error") {
    lam_optdesignsens <-
      optim(par = rep(0.01,K_asymp), fn = get_design_sens, psi_res = psi_YOverOmega_asymp, method = "Nelder-Mead", control = list(lmm = 5, pgtol = 1e-8, fnscale = -1e-3))$par
    }
  
  lamlist[[b]]=lam_optdesignsens
  
  dessensvec[b] <- get_design_sens(lamlist[[b]],psi_YOverOmega_asymp)
  
}

hist(dessensvec)
plot((unlist(lapply(lamlist,function(x)x[1]))),(unlist(lapply(lamlist,function(x)x[2]))))
coef(lm((unlist(lapply(lamlist,function(x)x[2]))) ~ (unlist(lapply(lamlist,function(x)x[1])))))
lamlist[[sort(dessensvec,index.return=T)$ix[nboot/2]]]





