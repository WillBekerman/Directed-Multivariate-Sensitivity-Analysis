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
nostratum = 175
# 
#tau_1 (impact on first outcome)
Tau_1 = -.5

#tau_2 (impact on second outcome)
Tau_2 = .25

#tau_3 (impact on third outcome)
Tau_3 = .5

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
  sv_whole_wholereused <- sv_analysis_analysisreused <- sv_analysis_randominit <- sv_analysis_randomnoise <- numeric(length=nsim)


for (sim in 1:nsim){
  #if (sim==6) browser()
  
  cat('\n\n\nSIMULATION NUMBER: ', sim, '\n\n\n')
  
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
  
  sensitivityResult_analysis_randominit=chiBarSquaredTest(Q = Q_analysis, #the data matrix
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
                                                               lam_init = runif(n=3,min=0,max=1))
  
  randomnoise = sensitivityResult_whole$OptLambda + rnorm(n=3,sd=.05)
  randomnoise[randomnoise < 0] = 0
  
  sensitivityResult_analysis_randomnoise=chiBarSquaredTest(Q = Q_analysis, #the data matrix
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
                                                          lam_init = randomnoise)
  
  
  res=list(
    sensitivityResult_whole=sensitivityResult_whole,
    sensitivityResult_planning=sensitivityResult_planning,
    sensitivityResult_analysis=sensitivityResult_analysis,
    sensitivityResult_analysis_withsearch=sensitivityResult_analysis_withsearch,
    sensitivityResult_analysis_analysisresused=sensitivityResult_analysis_analysisresused,
    sensitivityResult_whole_wholeresused=sensitivityResult_whole_wholeresused,
    sensitivityResult_analysis_randominit=sensitivityResult_analysis_randominit,
    sensitivityResult_analysis_randomnoise=sensitivityResult_analysis_randomnoise
  )
  
  if (any(unlist(lapply(res,function(i)i[[2]]))<1e-6)) {warning('There may be numerical problem \n')}
  
  sensvals <- lapply(res,function(i)i[[1]])
  sv_whole[sim] <- as.numeric(sensvals[1])
  sv_planning[sim] <- as.numeric(sensvals[2])
  sv_analysis[sim] <- as.numeric(sensvals[3])
  sv_analysis_withsearch[sim] <- as.numeric(sensvals[4])
  sv_analysis_analysisreused[sim] <- as.numeric(sensvals[5])
  sv_whole_wholereused[sim] <- as.numeric(sensvals[6])
  sv_analysis_randominit[sim] <- as.numeric(sensvals[7])
  sv_analysis_randomnoise[sim] <- as.numeric(sensvals[8])
}


gammas_vector <- seq(from=1,to=6,by=0.25)
nsim=78
power_sv_whole <- colSums(outer(sv_whole, gammas_vector, `>`))/nsim
power_sv_planning <- colSums(outer(sv_planning, gammas_vector, `>`))/nsim
power_sv_analysis <- colSums(outer(sv_analysis, gammas_vector, `>`))/nsim
power_sv_analysis_withsearch <- colSums(outer(sv_analysis_withsearch, gammas_vector, `>`))/nsim
power_sv_analysis_analysisreused <- colSums(outer(sv_analysis_analysisreused, gammas_vector, `>`))/nsim
power_sv_whole_wholereused <- colSums(outer(sv_whole_wholereused, gammas_vector, `>`))/nsim
power_sv_analysisrandominit <- colSums(outer(sv_analysis_randominit, gammas_vector, `>`))/nsim
power_sv_analysisrandomnoise <- colSums(outer(sv_analysis_randomnoise, gammas_vector, `>`))/nsim


methodnames <- c('Whole w/ Search','Analysis w/ Planning',
                 'Analysis w/ Search', 'Analysis w/ Analysis', 'Analysis w/ Random Init', 'Analysis w/ Analysis+Noise')
df <- data.frame(Gamma=rep(gammas_vector,length(methodnames)),
                 Sensitivity=c(power_sv_whole,
                               power_sv_analysis,
                               power_sv_analysis_withsearch,
                               power_sv_analysis_analysisreused,
                               power_sv_analysisrandominit,
                               power_sv_analysisrandomnoise),
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

#save(df,file='df_6dim_350sets.RData')
