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


Data_planning = syntheticData$normalData #this can be changed to use t-noised data
Data_planning$YData <- Data_planning$YData[ix_planning,]
Data_planning$s_k = apply(abs(Data_planning$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
psi_planning = psi_huber #this can be changed to use an inner-trimmed psi function

experimentalSetup_planning = makeExperimentalSetupFromYData(Data_planning, psi_planning)
Q_planning = experimentalSetup_planning$Q
populationSize_planning = dim(Q_planning)[2]
K_planning = dim(Q_planning)[1]
Z_planning = experimentalSetup_planning$Z #the treatment indicator
index_planning = experimentalSetup_planning$index

matchedSetAssignments_planning = rep(0, populationSize_planning)
for(ind in 1:length(index_planning))
{
  matchedSetAssignments_planning[unlist(index_planning[ind])] = ind
}


ptest_res = chiBarSquaredTest(Q = Q_planning, #the data matrix
                                               matchedSetAssignments = matchedSetAssignments_planning, #the stratum numbers for the individuals
                                               treatmentIndicator = Z_planning, #the treatment indicator
                                               numGamma = numGamma, #the number of Gammas to try
                                               alpha = .05, #the significance level of the test
                                               directions = directions, #the directions of each hypothesis
                                               step = 10, #optimization hyperparameter
                                               maxIter = 1000, #optimization hyperparameter
                                               showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                               verbose = verbose,
                                               outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Planning")





nboot=100
lamlist=vector('list',nboot)
dessensvec=numeric(nboot)
datmat=vector('list',nboot)
for (b in 1:nboot){
  Data_asymp = syntheticData$normalData
  Data_asymp$YData <- Data_asymp$YData[ix_planning,]
  
  bootix=sample(1:nrow(Data_asymp$YData),size=nrow(Data_asymp$YData),replace=T)
  bootdat=Data_asymp$YData[bootix,]
  Data_asymp$YData=bootdat
  Data_asymp$s_k = apply(abs(Data_asymp$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
  
  Data_planning=Data_asymp
  psi_planning = psi_huber #this can be changed to use an inner-trimmed psi function

  experimentalSetup_planning = makeExperimentalSetupFromYData(Data_planning, psi_planning)
  Q_planning = experimentalSetup_planning$Q
  populationSize_planning = dim(Q_planning)[2]
  K_planning = dim(Q_planning)[1]
  Z_planning = experimentalSetup_planning$Z #the treatment indicator
  index_planning = experimentalSetup_planning$index

  matchedSetAssignments_planning = rep(0, populationSize_planning)
  for(ind in 1:length(index_planning))
  {
    matchedSetAssignments_planning[unlist(index_planning[ind])] = ind
  }



  sensitivityResult_planning = chiBarSquaredTest(Q = Q_planning, #the data matrix
                                                 matchedSetAssignments = matchedSetAssignments_planning, #the stratum numbers for the individuals
                                                 treatmentIndicator = Z_planning, #the treatment indicator
                                                 numGamma = numGamma, #the number of Gammas to try
                                                 alpha = .05, #the significance level of the test
                                                 directions = directions, #the directions of each hypothesis
                                                 step = 10, #optimization hyperparameter
                                                 maxIter = 1000, #optimization hyperparameter
                                                 showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                 verbose = verbose,
                                                 outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Planning")


  lamlist[[b]]= sensitivityResult_planning$OptLambda #lam_optdesignsens

  dessensvec[b] <- sensitivityResult_planning$LargestRejectGamma#get_design_sens(lamlist[[b]],psi_YOverOmega_asymp)

  
  datmat[[b]] <- as.numeric(unlist(Data_planning$YData))

}

# hist(dessensvec)
# plot((unlist(lapply(lamlist,function(x)x[1]))),(unlist(lapply(lamlist,function(x)x[2]))))
# coef(lm((unlist(lapply(lamlist,function(x)x[2]))) ~ (unlist(lapply(lamlist,function(x)x[1])))))
# lamlist[[sort(dessensvec,index.return=T)$ix[nboot/2]]]


covmat=matrix(NA,3*nboot,nboot)
for(b in 1:nboot){
  covmat[,b] <- datmat[[b]]
}
covmat=cov(covmat)
invcovmat=solve(covmat)



#kap=6.5625/(1+6.5625) # this is planning smpl w/o bootstrapping w/ search
kap=ptest_res$LargestRejectGamma/(1+ptest_res$LargestRejectGamma)

#kap-sqrt(kap*(1-kap))

#kap-sqrt(kap*(1-kap))*sigq*((sqrt(qchibarsq(.95, diag(rep(1,3),3,3))))/(sqrt(.2))-(qnorm(.95))/((sqrt(.8))))-((dessensvec)/(1+dessensvec))
sigq=4/3#?
as.numeric(( kap-sqrt(kap*(1-kap))*sigq*((sqrt(qchibarsq(.95, diag(rep(1,3),3,3))))/(sqrt(.2))-(qnorm(.95))/((sqrt(.8))))-((dessensvec)/(1+dessensvec)) )^2 %*% invcovmat)







