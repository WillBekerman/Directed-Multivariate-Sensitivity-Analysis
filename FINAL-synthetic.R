#----- load libraries and source auxillary functions -----
setwd("C:/Users/bekerman/Documents/Github/Directed-Multivariate-Sensitivity-Analysis")
source(file = "./chiBarSquaredTest.R")
library(tidyverse)
#sources all files for generating synthetic data
allFiles = list.files(path = "./Synthetic Data Generation Functions")
for(file in allFiles)
{
  fileName = paste("./Synthetic Data Generation Functions/", file, sep = "")
  source(fileName)
}


#-----define simulation function -----
#Total number of Gamma's to try with the double-and-halve approach
#Epsilon (if the next Gamma differs from the current Gamma by < epsilon then quit)
#The significance level of the test
#Toggles printout amount
#Toggles diagnostics
#number of strata
# nostratum = 250
#impact on each outcome
# Taus <- rep(c(-0.5, -0.5, 0.25, 0.25, 0.5, 0.5),5)
#correlation
#directions
# directions=rep(c(rep("Less",2),rep("Greater",4)),5)
#split proportion
#number of simulations to average power

run_sim <- function(nostratum,Taus,directions,correlation,
                    other_baselines=c('Whole w/ Random',
                                      'Whole w/ Whole',
                                      'Whole w/ DS',
                                      'Analysis w/ Random',
                                      'Analysis w/ Analysis',
                                      'Analysis w/ DS'),
                    run_baseline=T, # Whole w/ Search
                    planning_sample_prop=0.2,nsim=100,seed=0,
                    numGamma=10,epsilon=1E-3,alpha=.05){
  
  # set seed for reproducibility
  set.seed(seed)
  
  #simulation metrics to return
  sv_whole <- sv_planning <- sv_analysis <- sv_whole_wholereused <- 
    sv_whole_randominit <- sv_whole_designsens <- sv_analysis_analysisreused <- 
    sv_analysis_randominit <- sv_analysis_designsens <- numeric(length=nsim)
  
  #split data into planning and analysis samples
  planning_sample_size_sets <- floor(planning_sample_prop*nostratum)
  ix_planning <- sort(sample(x=1:nostratum,size=planning_sample_size_sets))
  ix_analysis <- (1:nostratum)[-ix_planning]
  n1 <- planning_sample_size_sets
  n2 <- nostratum-planning_sample_size_sets
  
  #simulate design sensitivity from super-population
  syntheticData_asymp = generateData_general(rho = correlation, tauvec=Taus,nostratum = 1e7)
  syntheticData_asymp_small = generateData_general(rho = correlation, tauvec=Taus, nostratum = 1e5)
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
  #run if optimization didn't work
  if (class(lam_optdesignsens) == "try-error") {
    lam_optdesignsens <-
      optim(par = rep(0.01,K_asymp), fn = get_design_sens, psi_res = psi_YOverOmega_asymp, method = "L-BFGS-B", lower = 0, upper = 1, control = list(lmm = 5, pgtol = 1e-8, fnscale = -1))$par
  }
  
  
  #-----loop over simulations-----
  for (sim in 1:nsim){
    
    cat('\n\n\n\n\nSIMULATION NUMBER: ', sim, '\n\n\n\n\n\n')
    
    #-----create and partition synthetic data-----
    syntheticData = generateData_general(rho = correlation, tauvec=Taus, nostratum = nostratum)
    
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
    
    
    #-----run sensitivity analyses-----
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
    
    ## BASELINE WHOLE W/ SEARCH ##
    if (run_baseline){
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
    }

    ## OTHER BASELINES ##
    if ('Whole w/ Random' %in% other_baselines){
      randominit = runif(n=length(Taus))
      sensitivityResult_whole_randominit=chiBarSquaredTest(Q = Q_whole, #the data matrix
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
                                                           lam_init = randominit)
    }
    if ('Whole w/ Whole' %in% other_baselines){
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
    }
    if ('Whole w/ DS' %in% other_baselines){
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
    }
    
    if ('Analysis w/ Random' %in% other_baselines){
      randominit = runif(n=length(Taus))
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
                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Analysis",
                                                              lam_init = randominit)
    }
    if ('Analysis w/ Analysis' %in% other_baselines){
      sensitivityResult_analysis_withsearch=chiBarSquaredTest(Q = Q_analysis, #the data matrix
                                                              matchedSetAssignments = matchedSetAssignments_analysis, #the stratum numbers for the individuals
                                                              treatmentIndicator = Z_analysis, #the treatment indicator
                                                              numGamma = numGamma, #the number of Gammas to try
                                                              alpha = .05, #the significance level of the test
                                                              directions = directions, #the directions of each hypothesis
                                                              step = 10, #optimization hyperparameter
                                                              maxIter = 1000, #optimization hyperparameter
                                                              showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                                              verbose = verbose,
                                                              outputDirName = "SyntheticExample_Sensitivity_Analysis_Results_Wnole")
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
    }
    if ('Analysis w/ DS' %in% other_baselines){
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
      
    
    }
    
    
    #----update metrics-----
    sv_planning[sim] <- sensitivityResult_planning$LargestRejectGamma
    sv_analysis[sim] <- sensitivityResult_analysis$LargestRejectGamma
     
    if(run_baseline){
      sv_whole[sim] <- sensitivityResult_whole$LargestRejectGamma
    }
    if ('Whole w/ Random' %in% other_baselines){
      sv_whole_randominit[sim] <- 
        sensitivityResult_whole_randominit$LargestRejectGamma
    }
    if ('Whole w/ Whole' %in% other_baselines){
      sv_whole_wholereused[sim] <- 
        sensitivityResult_whole_wholeresused$LargestRejectGamma
    }
    if ('Whole w/ DS' %in% other_baselines){
      sv_whole_designsens[sim] <- 
        sensitivityResult_whole_designsens$LargestRejectGamma
    }
    if ('Analysis w/ Random' %in% other_baselines){
      sv_analysis_randominit[sim] <- 
        sensitivityResult_analysis_randominit$LargestRejectGamma
    }
    if ('Analysis w/ Analysis' %in% other_baselines){
      sv_analysis_analysisreused[sim] <- 
        sensitivityResult_analysis_analysisresused$LargestRejectGamma
    }
    if ('Analysis w/ DS' %in% other_baselines){
      sv_analysis_designsens[sim] <- 
        sensitivityResult_analysis_designsens$LargestRejectGamma
    }

  }#end of sim
  
  #return metrics as a list called res
  res=list()
  #res$sv_planning=sv_planning
  res$sv_analysis=sv_analysis
  res$sv_whole=sv_whole
  res$sv_analysis_analysisreused=sv_analysis_analysisreused
  res$sv_analysis_designsens=sv_analysis_designsens
  res$sv_analysis_randominit=sv_analysis_randominit
  res$sv_whole_wholereused=sv_whole_wholereused
  res$sv_whole_designsens=sv_whole_designsens
  res$sv_whole_randominit=sv_whole_randominit
  #remove only zeros (NA entries)
  res=res[unlist(lapply(res,function(x)!all(x==0)))]

  return(res)
  
}


#-----run simulations -----
#set desired parameters
# verbose=T
# showDiagnostics=T
# nsim=3
# nostratum=200
# Taus=c(-.5,.25,.5)
# directions=c('Less','Greater','Greater')
# corr=0
# other_baselines <- c('Whole w/ DS')#,'Analysis w/ DS')
# methodnames <- c('Analysis w/ Planning','Whole w/ Search', other_baselines)

verbose=T
showDiagnostics=T
nsim=1000
nostratum=1000
Taus=rep(c(-.5,.25,.5,.75,-.75),7)
directions=rep(c('Less','Greater','Greater','Greater','Less'),7)
corr=0
other_baselines <- c('Whole w/ DS')#,'Analysis w/ DS')
methodnames <- c('Analysis w/ Planning','Whole w/ Search', other_baselines)


sim_result=run_sim(nostratum=nostratum,Taus=Taus,directions=directions,correlation=corr,
                   other_baselines=other_baselines,run_baseline=T,
                   planning_sample_prop=0.2,
                   nsim=nsim)


# gammas_vector <- seq(from=1,to=5,by=0.25)
# gammas_vectorold <- seq(from=1,to=5.6,by=.92/4)
gammas_vectorold <- seq(from=1,to=6,by=.01)
gammas_vector <- exp(gammas_vectorold)
power_result <- lapply(sim_result,function(x)colSums(outer(x,gammas_vector,`>`))/nsim)

df <- data.frame(Gamma=rep(gammas_vectorold,length(methodnames)),
                 Sensitivity=unlist(power_result),
                 Method=rep(methodnames,each=length(gammas_vector)))


# ggplot(df, aes(x = Gamma, y = Sensitivity)) +
#   geom_line(aes(col=Method), size=1) +
#   labs(
#     x = expression(gamma),#str2expression(paste0('log(',expression(Gamma),')')),
#     y = "Average Power"
#   ) + theme_classic() + theme(legend.position="none")
# 








# Create a function to smooth and clip the data
smooth_and_clip <- function(df, xvar, yvar, method = "scam", span = 0.5) {
  # Apply the smoothing function
  smoothed <- df %>%
    group_by(Method) %>%
    mutate(Smoothed = if (method == "loess") {
    predict(loess(get(yvar) ~ get(xvar), span = span))
  } else if (method == "lm") {
    predict(lm(get(yvar) ~ get(xvar)))
  } else if (method == "scam") {
    predict( scam(get(yvar) ~ s(get(xvar), k = 25, bs = "mpd")) )
  })
  # Clip the smoothed values between 0 and 1
  smoothed <- smoothed %>%
    mutate(Smoothed = pmax(pmin(Smoothed, 1), 0))
  
  return(smoothed)
}

# Apply the function to your data
df_smoothed <- smooth_and_clip(df, xvar = "Gamma", yvar = "Sensitivity", method = "scam")



ggplot(df_smoothed, aes(x = Gamma, y = Smoothed, col=Method)) +
  geom_line(size = 1.25) +
  ylim(-0.001, 1.001)+#xlim(1,6)+
  labs(
    x = expression(gamma),#str2expression(paste0('log(',expression(Gamma),')')),
    y = "Average Power"
  ) + theme_classic() + theme(legend.position="none") + lims(x=c(1,5.8))



# ggplot(df, aes(x = Gamma, y = Sensitivity, col=Method)) +
#   geom_line(size = 1.25) +
#   ylim(-0.01, 1.01)+#xlim(1,6)+
#   labs(
#     x = expression(gamma),#str2expression(paste0('log(',expression(Gamma),')')),
#     y = "Average Power"
#   ) + theme_classic() + theme(legend.position="none")

library(rlang)
df_smoothed$Method[df_smoothed$Method=="Whole w/ Search"]=paste0('Whole w/ Scheffe')
# df_smoothed$Method[df_smoothed$Method=="Whole w/ DS"]=exprs(expression(paste0('Whole w/ ', (lambda),'DS')))
# df_smoothed$Method[df_smoothed$Method=="Analysis w/ Planning"]=exprs(expression(paste0('Analysis w/ ', (lambda),'Planning')))



# df_smoothed <- df_smoothed %>% mutate(
#   Method = (case_when(
#     Method=="Whole w/ Search" ~ 'Whole w/ Scheffe',
#     Method=="Whole w/ DS" ~ (exprs(paste0('Whole w/ ', (lambda),'DS'))),
#     Method=="Analysis w/ Planning" ~ (exprs(paste0('Analysis w/ ', (lambda),'Planning'))),
#   ))
# )


  
library(cowplot)
get_legend_plot <- ggplot(df_smoothed, aes(x = Gamma, y = Smoothed, col=Method)) +
  geom_line(size = 1.25) +
  ylim(-0.001, 1.001)+#xlim(1,6)+
  labs(
    x = expression(gamma),#str2expression(paste0('log(',expression(Gamma),')')),
    y = "Average Power"
  ) + theme_classic() + lims(x=c(1,5.8))+theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

all_plots_nuc <- simrerssmall_plt
all_plots_uc <- simrersbig_plt

legend<- get_plot_component(get_legend_plot, 'guide-box-bottom', return_all = TRUE)
combined_plot <- plot_grid(NULL, all_plots_nuc, NULL, all_plots_uc, legend, ncol=1, rel_heights = c(.1,1,.1,1,.1), labels = c('', 'I = 250', '', 'I = 1000', ''), vjust=-0.5)

combined_plot



