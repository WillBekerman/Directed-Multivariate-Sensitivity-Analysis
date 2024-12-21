################################################################################
#           An example of how to use the chibarsq code on synthetic data
#           Follows the simulation set-up from Cohen, Olson, and Fogart (2018)
#                 Generates trivariate data at Gamma = 1 and performs sensitivity analysis on it
################################################################################
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
nostratum = 100
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

#number of simulations to average power
nsim=100

#sim metrics to return
sv_whole <- sv_planning <- sv_analysis <- sv_analysis_withsearch <- sv_whole_planreused <- 
  sv_whole_wholereused <- sv_analysis_analysisreused <- sv_analysis_randominit <- sv_analysis_randomnoise <- 
  sv_analysis_randomnoise_.01 <- sv_analysis_randomnoise_.05 <- sv_analysis_randomnoise_.1 <- sv_analysis_randomnoise_.2 <- 
  sv_whole_designsens <- sv_analysis_designsens <- numeric(length=nsim)


syntheticData_asymp = generateData3(rho = correlation, Tau_1 = Tau_1, Tau_2 = Tau_2, Tau_3 = Tau_3, nostratum = 1e7)
syntheticData_asymp_small = generateData3(rho = correlation, Tau_1 = Tau_1, Tau_2 = Tau_2, Tau_3 = Tau_3, nostratum = 1e6)
omega_k <- apply(abs(syntheticData_asymp$normalData$YData), MARGIN = 2, FUN = median)

Data_asymp = syntheticData_asymp_small$normalData
Data_asymp$s_k = apply(abs(Data_asymp$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
Data_asymp$omega_k = omega_k  # omega_k is set to the median currently
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
  lam_optdesignsens <-
    optim(par = rep(0.01,K_asymp), fn = get_design_sens, psi_res = psi_YOverOmega_asymp, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1))$par
}



for (sim in 1:nsim){

  cat('\n\n\n\n\nSIMULATION NUMBER: ', sim, '\n\n\n\n\n\n')
  
  ################################################################################
  #                         Creates the synthetic data
  ################################################################################
  syntheticData = generateData3(rho = correlation, Tau_1 = Tau_1, Tau_2 = Tau_2, Tau_3 = Tau_3, nostratum = nostratum)
  #syntheticData = generateData_general(rho = correlation, tauvec=Taus, nostratum = nostratum)
  
  ################################################################################
  #                           We take in the data and turn it into
  #                           a format that works for the chibarsq-test code.
  ################################################################################
  
  ### WHOLE SAMPLE ###
  Data_whole = syntheticData$normalData #this can be changed to use t-noised data
  Data_whole$s_k = apply(abs(Data_whole$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
  Data_whole$omega_k = omega_k  # omega_k is set to the median currently
  psi_whole = psi_huber #this can be changed to use an inner-trimmed psi function
  
  experimentalSetup_whole = makeExperimentalSetupFromYData(Data_whole, psi_whole)
  Q_whole = experimentalSetup_whole$Q
  populationSize_whole = dim(Q_whole)[2]
  K_whole = dim(Q_whole)[1]
  Z_whole = experimentalSetup_whole$Z #the treatment indicator
  index_whole = experimentalSetup_whole$index
  
  matchedSetAssignments_whole = rep(0, populationSize_whole)
  for(ind in 1:length(index_whole))
  {
    matchedSetAssignments_whole[unlist(index_whole[ind])] = ind
  }
  
  ### PLANNING SAMPLE ###
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
  
  ### ANALYSIS SAMPLE ###
  Data_analysis = syntheticData$normalData #this can be changed to use t-noised data
  Data_analysis$YData <- Data_analysis$YData[ix_analysis,]
  Data_analysis$s_k = apply(abs(Data_analysis$YData), MARGIN = 2, FUN = median)  # s_k is set to the median currently
  #Data_analysis$omega_k = apply(abs(syntheticData_asymp$normalData$YData), MARGIN = 2, FUN = median)  # omega_k is set to the median currently
  psi_analysis = psi_huber #this can be changed to use an inner-trimmed psi function
  
  experimentalSetup_analysis = makeExperimentalSetupFromYData(Data_analysis, psi_analysis)
  Q_analysis = experimentalSetup_analysis$Q
  populationSize_analysis = dim(Q_analysis)[2]
  K_analysis = dim(Q_analysis)[1]
  Z_analysis = experimentalSetup_analysis$Z #the treatment indicator
  index_analysis = experimentalSetup_analysis$index
  
  matchedSetAssignments_analysis = rep(0, populationSize_analysis)
  for(ind in 1:length(index_analysis))
  {
    matchedSetAssignments_analysis[unlist(index_analysis[ind])] = ind
  }
  
  # # D = syntheticData_asymp_small$normalData
  # # D$omega_k = Data_analysis$omega_k
  # # experimentalSetup_analysis_asymp = makeExperimentalSetupFromYData_asymp(D, psi_analysis)
  # # psi_YOverOmega_analysis = experimentalSetup_analysis_asymp$psi_YOverOmega
  # # psi_YOverOmega_analysis[,directions == 'Less']=-psi_YOverOmega_analysis[,directions == 'Less']
  # # # experimentalSetup_analysis_asymp = makeExperimentalSetupFromYData_asymp(Data_analysis, psi_analysis)
  # # # psi_YOverOmega_analysis = experimentalSetup_analysis_asymp$psi_YOverOmega
  # # # psi_YOverOmega_analysis[,directions == 'Less']=-psi_YOverOmega_analysis[,directions == 'Less']
  # # eta <- function(lam,psi_res){
  # #   # mean(abs( lam%*%t(psi_res) ))
  # #   mean(sqrt( (lam%*%t(psi_res))^2 ))
  # # }
  # # theta <- function(lam,psi_res){
  # #   mean( lam%*%t(psi_res) )
  # # }
  # # get_design_sens <- function(lam,psi_res){
  # #   (eta(lam,psi_res)+theta(lam,psi_res)) / (eta(lam,psi_res)-theta(lam,psi_res))
  # # }
  # # #lam_optdesignsens_analysis = optim(par = rep(0.01,K_analysis), fn = get_design_sens, psi_res = psi_YOverOmega_analysis, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1e-3))$par
  # # lam_optdesignsens_analysis <- try({
  # #   optim(par = rep(0.01,K_analysis), fn = get_design_sens, psi_res = psi_YOverOmega_analysis, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1e-3))$par
  # # })
  # # if (class(lam_optdesignsens_analysis) == "try-error") {
  # #   lam_optdesignsens_analysis <- 
  # #     optim(par = rep(0.01,K_analysis), fn = get_design_sens, psi_res = psi_YOverOmega_analysis, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1))$par
  # # }
  # # cat('\n\nDESIGN SENSITIVITY OPT LAM:', lam_optdesignsens_analysis, '\n\n')
  # lam_optdesignsens_analysis=c(0.02,0.12,0.02) ###### THIS IS WRONG!!!!!!
  
  ################################################################################
  #                   Perform the Sensitivity Analysis
  ################################################################################
  
  ### WHOLE SAMPLE ###
  sensitivityResult_whole = chiBarSquaredTest(Q = Q_whole, #the data matrix
                                              matchedSetAssignments = matchedSetAssignments_whole, #the stratum numbers for the individuals
                                              treatmentIndicator = Z_whole, #the treatment indicator
                                              numGamma = numGamma, #the number of Gammas to try
                                              alpha = .05, #the significance level of the test
                                              directions = directions, #the directions of each hypothesis
                                              step = 10, #optimization hyperparameter
                                              maxIter = 1000, #optimization hyperparameter
                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                              verbose = verbose,
                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Whole")
  
  
  ### PLANNING SAMPLE ###
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
  
  ### ANALYSIS SAMPLE ###
  sensitivityResult_analysis = chiBarSquaredTest(Q = Q_analysis, #the data matrix
                                                 matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
                                                 treatmentIndicator = Z_analysis, #the treatment indicator
                                                 numGamma = numGamma, #the number of Gammas to try
                                                 alpha = .05, #the significance level of the test
                                                 directions = directions, #the directions of each hypothesis
                                                 step = 10, #optimization hyperparameter
                                                 maxIter = 1000, #optimization hyperparameter
                                                 showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                 verbose = verbose,
                                                 outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Analysis",
                                                 lam_init=sensitivityResult_planning$OptLambda)
  
  
  ## SOME COMPARISONS ##
  
  sensitivityResult_analysis_withsearch = chiBarSquaredTest(Q = Q_analysis, #the data matrix
                                                            matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
                                                            treatmentIndicator = Z_analysis, #the treatment indicator
                                                            numGamma = numGamma, #the number of Gammas to try
                                                            alpha = .05, #the significance level of the test
                                                            directions = directions, #the directions of each hypothesis
                                                            step = 10, #optimization hyperparameter
                                                            maxIter = 1000, #optimization hyperparameter
                                                            showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                            verbose = verbose,
                                                            outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Analysis")
  
  
  sensitivityResult_analysis_analysisresused=chiBarSquaredTest(Q = Q_analysis, #the data matrix
                                                               matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
                                                               treatmentIndicator = Z_analysis, #the treatment indicator
                                                               numGamma = numGamma, #the number of Gammas to try
                                                               alpha = .05, #the significance level of the test
                                                               directions = directions, #the directions of each hypothesis
                                                               step = 10, #optimization hyperparameter
                                                               maxIter = 1000, #optimization hyperparameter
                                                               showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                               verbose = verbose,
                                                               outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
                                                               lam_init = sensitivityResult_analysis_withsearch$OptLambda)
  
  
  sensitivityResult_whole_wholeresused=chiBarSquaredTest(Q = Q_whole, #the data matrix
                                                         matchedSetAssignments = matchedSetAssignments_whole, #the stratum numbers for the individuals
                                                         treatmentIndicator = Z_whole, #the treatment indicator
                                                         numGamma = numGamma, #the number of Gammas to try
                                                         alpha = .05, #the significance level of the test
                                                         directions = directions, #the directions of each hypothesis
                                                         step = 10, #optimization hyperparameter
                                                         maxIter = 1000, #optimization hyperparameter
                                                         showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                         verbose = verbose,
                                                         outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
                                                         lam_init = sensitivityResult_whole$OptLambda)
  # 
  # randomnoise = sensitivityResult_whole$OptLambda + rnorm(n=3,sd=.01)
  # randomnoise[randomnoise < 0] = 0
  # sensitivityResult_analysis_randomnoise_.01=chiBarSquaredTest(Q = Q_analysis, #the data matrix
  #                                                          matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
  #                                                          treatmentIndicator = Z_analysis, #the treatment indicator
  #                                                          numGamma = numGamma, #the number of Gammas to try
  #                                                          alpha = .05, #the significance level of the test
  #                                                          directions = directions, #the directions of each hypothesis
  #                                                          step = 10, #optimization hyperparameter
  #                                                          maxIter = 1000, #optimization hyperparameter
  #                                                          showDiagnostics = showDiagnostics, #whether or not diagonstics are output
  #                                                          verbose = verbose,
  #                                                          outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
  #                                                          lam_init = randomnoise)
  # 
  # 
  # 
  # 
  # randomnoise = sensitivityResult_whole$OptLambda + rnorm(n=3,sd=.05)
  # randomnoise[randomnoise < 0] = 0
  # sensitivityResult_analysis_randomnoise_.05=chiBarSquaredTest(Q = Q_analysis, #the data matrix
  #                                                              matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
  #                                                              treatmentIndicator = Z_analysis, #the treatment indicator
  #                                                              numGamma = numGamma, #the number of Gammas to try
  #                                                              alpha = .05, #the significance level of the test
  #                                                              directions = directions, #the directions of each hypothesis
  #                                                              step = 10, #optimization hyperparameter
  #                                                              maxIter = 1000, #optimization hyperparameter
  #                                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
  #                                                              verbose = verbose,
  #                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
  #                                                              lam_init = randomnoise)
  # 
  # 
  # 
  # randomnoise = sensitivityResult_whole$OptLambda + rnorm(n=3,sd=.1)
  # randomnoise[randomnoise < 0] = 0
  # sensitivityResult_analysis_randomnoise_.1=chiBarSquaredTest(Q = Q_analysis, #the data matrix
  #                                                              matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
  #                                                              treatmentIndicator = Z_analysis, #the treatment indicator
  #                                                              numGamma = numGamma, #the number of Gammas to try
  #                                                              alpha = .05, #the significance level of the test
  #                                                              directions = directions, #the directions of each hypothesis
  #                                                              step = 10, #optimization hyperparameter
  #                                                              maxIter = 1000, #optimization hyperparameter
  #                                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
  #                                                              verbose = verbose,
  #                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
  #                                                              lam_init = randomnoise)
  # 
  # 
  # 
  # randomnoise = sensitivityResult_whole$OptLambda + rnorm(n=3,sd=.2)
  # randomnoise[randomnoise < 0] = 0
  # sensitivityResult_analysis_randomnoise_.2=chiBarSquaredTest(Q = Q_analysis, #the data matrix
  #                                                              matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
  #                                                              treatmentIndicator = Z_analysis, #the treatment indicator
  #                                                              numGamma = numGamma, #the number of Gammas to try
  #                                                              alpha = .05, #the significance level of the test
  #                                                              directions = directions, #the directions of each hypothesis
  #                                                              step = 10, #optimization hyperparameter
  #                                                              maxIter = 1000, #optimization hyperparameter
  #                                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
  #                                                              verbose = verbose,
  #                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
  #                                                              lam_init = randomnoise)
  
  
  
  sensitivityResult_whole_designsens=chiBarSquaredTest(Q = Q_whole, #the data matrix
                                                              matchedSetAssignments = matchedSetAssignments_whole, #the stratum numbers for the individuals
                                                              treatmentIndicator = Z_whole, #the treatment indicator
                                                              numGamma = numGamma, #the number of Gammas to try
                                                              alpha = .05, #the significance level of the test
                                                              directions = directions, #the directions of each hypothesis
                                                              step = 10, #optimization hyperparameter
                                                              maxIter = 1000, #optimization hyperparameter
                                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                              verbose = verbose,
                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
                                                              lam_init = lam_optdesignsens)
  
  sensitivityResult_analysis_designsens=chiBarSquaredTest(Q = Q_analysis, #the data matrix
                                                       matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
                                                       treatmentIndicator = Z_analysis, #the treatment indicator
                                                       numGamma = numGamma, #the number of Gammas to try
                                                       alpha = .05, #the significance level of the test
                                                       directions = directions, #the directions of each hypothesis
                                                       step = 10, #optimization hyperparameter
                                                       maxIter = 1000, #optimization hyperparameter
                                                       showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                       verbose = verbose,
                                                       outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole",
                                                       lam_init = lam_optdesignsens)
  
  
  
  
  res=list(
    # sensitivityResult_whole=sensitivityResult_whole,
    sensitivityResult_planning=sensitivityResult_planning,
    sensitivityResult_analysis=sensitivityResult_analysis,
    sensitivityResult_analysis_analysisresused=sensitivityResult_analysis_analysisresused,
    sensitivityResult_whole_wholeresused=sensitivityResult_whole_wholeresused,
    # sensitivityResult_analysis_randomnoise_.01=sensitivityResult_analysis_randomnoise_.01,
    # sensitivityResult_analysis_randomnoise_.05=sensitivityResult_analysis_randomnoise_.05,
    # sensitivityResult_analysis_randomnoise_.1=sensitivityResult_analysis_randomnoise_.1,
    # sensitivityResult_analysis_randomnoise_.2=sensitivityResult_analysis_randomnoise_.2,
    sensitivityResult_whole_designsens=sensitivityResult_whole_designsens,
    sensitivityResult_analysis_designsens=sensitivityResult_analysis_designsens
  )
  
  #if (any(unlist(lapply(res,function(i)i[[2]]))<1e-6)) {warning('There may be numerical problem \n')}
  
  sensvals <- lapply(res,function(i)i[[1]])
  # sv_whole[sim] <- as.numeric(sensvals[1])
  sv_planning[sim] <- as.numeric(sensvals[1])
  sv_analysis[sim] <- as.numeric(sensvals[2])
  sv_analysis_analysisreused[sim] <- as.numeric(sensvals[3])
  sv_whole_wholereused[sim] <- as.numeric(sensvals[4])
  # sv_analysis_randomnoise_.01[sim] <- as.numeric(sensvals[6])
  # sv_analysis_randomnoise_.05[sim] <- as.numeric(sensvals[7])
  # sv_analysis_randomnoise_.1[sim] <- as.numeric(sensvals[8])
  # sv_analysis_randomnoise_.2[sim] <- as.numeric(sensvals[9])
  sv_whole_designsens[sim] <- as.numeric(sensvals[5])
  sv_analysis_designsens[sim] <- as.numeric(sensvals[6])
}


gammas_vector <- seq(from=1,to=10,by=0.25)
# power_sv_whole <- colSums(outer(sv_whole, gammas_vector, `>`))/nsim
power_sv_planning <- colSums(outer(sv_planning, gammas_vector, `>`))/nsim
power_sv_analysis <- colSums(outer(sv_analysis, gammas_vector, `>`))/nsim
power_sv_analysis_analysisreused <- colSums(outer(sv_analysis_analysisreused, gammas_vector, `>`))/nsim
power_sv_whole_wholereused <- colSums(outer(sv_whole_wholereused, gammas_vector, `>`))/nsim
# power_sv_analysisrandomnoise_.01 <- colSums(outer(sv_analysis_randomnoise_.01, gammas_vector, `>`))/nsim
# power_sv_analysisrandomnoise_.05 <- colSums(outer(sv_analysis_randomnoise_.05, gammas_vector, `>`))/nsim
# power_sv_analysisrandomnoise_.1 <- colSums(outer(sv_analysis_randomnoise_.1, gammas_vector, `>`))/nsim
# power_sv_analysisrandomnoise_.2 <- colSums(outer(sv_analysis_randomnoise_.2, gammas_vector, `>`))/nsim
power_sv_whole_designsens <- colSums(outer(sv_whole_designsens, gammas_vector, `>`))/nsim
power_sv_analysis_designsens <- colSums(outer(sv_analysis_designsens, gammas_vector, `>`))/nsim


# 
# methodnames <- c('Whole w/ Search','Analysis w/ Planning',
#                  'Analysis w/ Analysis', 'Analysis w/ Analysis+Noise(0.01)',
#                  'Analysis w/ Analysis+Noise(0.05)','Analysis w/ Analysis+Noise(0.1)','Analysis w/ Analysis+Noise(0.2)',
#                  'Analysis w/ Max DesSens')
# df <- data.frame(Gamma=rep(gammas_vector,length(methodnames)),
#                  Sensitivity=c(power_sv_whole,
#                                power_sv_analysis,
#                                power_sv_analysis_analysisreused,
#                                power_sv_analysisrandomnoise_.01,
#                                power_sv_analysisrandomnoise_.05,
#                                power_sv_analysisrandomnoise_.1,
#                                power_sv_analysisrandomnoise_.2,
#                                power_sv_analysis_designsens),
#                  Method=rep(methodnames,each=length(gammas_vector)))


methodnames <- c('Analysis w/ Planning','Whole w/ Whole',
                 'Analysis w/ Analysis', 'Whole w/ Max DesSens', 'Analysis w/ Max DesSens')
df <- data.frame(Gamma=rep(gammas_vector,length(methodnames)),
                 Sensitivity=c(power_sv_analysis,
                               power_sv_whole_wholereused,
                               power_sv_analysis_analysisreused,
                               power_sv_whole_designsens,
                               power_sv_analysis_designsens),
                 Method=rep(methodnames,each=length(gammas_vector)))


ggplot(df, aes(x = Gamma, y = Sensitivity)) +
  geom_line(aes(linetype=Method,col=Method), size=1) +
  labs(
    x = expression(Gamma),
    y = "Average Power"
  )
ggplot(df, aes(x = Gamma, y = Sensitivity)) +
  geom_smooth(aes(linetype=Method,col=Method), size=1,se = FALSE) + 
  labs(
    x = expression(Gamma),
    y = "Average Power"
  )

#save(df,file='df_noisy_3dim_175sets.RData')
