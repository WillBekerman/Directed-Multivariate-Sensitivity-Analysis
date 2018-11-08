################################################################################
#           An example of how to use the chibarsq code on synthetic data
#           Follows the simulation set-up from Cohen, Olson, and Fogart (2018)
#                 Generates trivariate data at Gamma = 1 and performs sensitivity analysis on it
################################################################################

source(file = "./chiBarSquaredTest.R")

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
verbose = TRUE

#Toggles diagnostics
showDiagnostics = TRUE

#number of strata
nostratum = 100

#tau_1 (impact on first outcome)
Tau_1 = .5

#tau_2 (impact on second outcome)
Tau_2 = -.5

#tau_3 (impact on third outcome)
Tau_3 = -.5

#correlation
correlation = 0

#directions
directions = c("Greater", "Less", "Less")
################################################################################
#                         Creates the synthetic data
################################################################################
syntheticData = generateData3(rho = correlation, Tau_1 = Tau_1, Tau_2 = Tau_2, Tau_3 = Tau_3, nostratum = nostratum)

Data = syntheticData$normalData #this can be changed to use t-noised data
psi = psi_huber #this can be changed to use an inner-trimmed psi function

################################################################################
#                           We take in the data and turn it into
#                           a format that works for the chibarsq-test code.
################################################################################
experimentalSetup = makeExperimentalSetupFromYData(Data, psi)
Q = experimentalSetup$Q
populationSize = dim(Q)[2]
K = dim(Q)[1]
Z = experimentalSetup$Z #the treatment indicator
index = experimentalSetup$index


matchedSetAssignments = rep(0, populationSize)
for(ind in 1:length(index))
{
  matchedSetAssignments[unlist(index[ind])] = ind
}

################################################################################
#                   Perform the Sensitivity Analysis
################################################################################
sensitivityResult = chiBarSquaredTest(Q = Q, #the data matrix
                                      matchedSetAssignments = matchedSetAssignments, #the stratum numbers for the individuals
                                      treatmentIndicator = Z, #the treatment indicator
                                      numGamma = numGamma, #the number of Gammas to try
                                      alpha = .05, #the significance level of the test
                                      directions = directions, #the directions of each hypothesis
                                      step = 100, #optimization hyperparameter
                                      maxIter = 1000, #optimization hyperparameter
                                      showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                      verbose = verbose,
                                      outputDirName = "SyntheticExample_Sensitivity_Analysis_Results")