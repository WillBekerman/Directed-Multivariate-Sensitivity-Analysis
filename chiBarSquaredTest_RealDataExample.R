################################################################################
#               An example of how to use the chibarsq code on real data
################################################################################

################################################################################
#                                   Set Up
################################################################################
# install.packages('DOS') #Install the 'DOS' package (if you have not already done) so by uncommenting an running this line
library('DOS') #The R package to "Design of Observational Studies" (2010) by P. Rosenbaum
source(file = "./chiBarSquaredTest.R")

################################################################################
#                     The example that we will use is data from:
#                     Angrist, J. D. and Lavy, V. (1999). Using Maimonides' rule to estimate the effect of class size on scholastic achievement. The Quarterly Journal of Economics, 114, 533-575.

#                     To learn more about this data-set enter "?angristlavy" into the console.
################################################################################

data("angristlavy")

################################################################################
#                             Tweakable Parameters
################################################################################
#Total number of Gamma's to try with the double-and-halve approach
numGamma = 15

#Epsilon (if the next Gamma differs from the current Gamma by < epsilon then quit)
epsilon = 1E-3

#The significance level of the test
alpha = .05

#Toggles printout amount
verbose = TRUE

#Toggles diagnostics
showDiagnostics = FALSE


################################################################################
#                     Extracts Useful information
################################################################################
matchedSetAssignments = angristlavy$pair #the stratum that each person was assigned to
Z = 1 - angristlavy$z #treatment indicator

rawData = data.frame(math = angristlavy$avgmath, verb = angristlavy$avgverb)

index = list() #the ith element of the list is the raw indices of the individuals in stratum i.
for(ind in 1:length(unique(matchedSetAssignments)))
{
  index[[ind]] = which(matchedSetAssignments == ind)
}

################################################################################
#                         Processes the data:
#                         We will use a Wilcoxon signed rank statistic to form the Q matrix
################################################################################
Q = matrix(0, nrow = dim(rawData)[1], ncol = dim(rawData)[2])
treatmentMinusContol = matrix(0, nrow = (dim(rawData)[1]/2), ncol = dim(rawData)[2])

for(strat in unique(matchedSetAssignments)) #iterates over strata in the data-set
{
  treatmentIndex = base::which(Z[unlist(index[strat])] == 1)
  controlIndex = base::which(Z[unlist(index[strat])] == 0)
  treatmentMinusContol[strat, ] = t(rawData[unlist(index[strat])[treatmentIndex], ] - rawData[unlist(index[strat])[controlIndex], ])
}

rankedTreatmentMinusControl = apply(abs(treatmentMinusContol), MARGIN = 2, rank)

for(individual in 1:dim(rawData)[1]) #iterates over individuals in the data-set
{
  strat = matchedSetAssignments[individual]
  indivsData = rawData[individual, ]
  pairsData = rawData[unlist(index[strat])[unlist(index[strat]) != individual], ]
  
  s = as.numeric(indivsData > pairsData)
  
  Q[individual,] = s * rankedTreatmentMinusControl[strat, ]
}

Q = t(Q) #We need to transpose the Q matrix to fit the format of the chiBarSquaredTest() function

################################################################################
#                         Perform the Sensitivity Analysis (with a greater-than alternative on each outcome)
################################################################################
sensitivityResult = chiBarSquaredTest(Q = Q, #the data matrix
                                matchedSetAssignments = matchedSetAssignments, #the stratum numbers for the individuals
                                treatmentIndicator = Z, #the treatment indicator
                                numGamma = 15, #the number of Gammas to try
                                alpha = .05, #the significance level of the test
                                directions = "Greater", #the directions of each hypothesis
                                step = 100, #optimization hyperparameter
                                maxIter = 1000, #optimization hyperparameter
                                showDiagnostics = showDiagnostics, #whether or not diagonstics are output
                                verbose = verbose,
                                outputDirName = "Example_Sensitivity_Analysis_Results")