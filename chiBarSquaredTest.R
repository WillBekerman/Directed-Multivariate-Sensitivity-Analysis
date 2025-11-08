################################################################################
#                         Chi-Bar-Squared Test 
#                         (Uses the method of Cohen, Olson, and Fogarty 2018)
################################################################################
#' Chi-Bar-Squared Test
#' 
#' Implementation of the method of Cohen, Olson, and Fogarty 2018
#' 
#' @param Q the matrix of the q_{ijk}.  It has K rows (the number of outcomes) and N columns (the number of individuals)
#' @param matchedSetAssignments the vector where the t^th element is the is the number of the matched set for the t^th individual
#' @param treatmentIndicator the vector where the t^th element is the indicator of treatment for the t^th individual (0 if control, 1 if treated)
#' @param numGamma the number of Gamma to try in a double-and-halve sensitivity analysis (defaults to 10)
#' @param alpha the significance level of the test (defaults to .05)
#' @param directions the vector where the k^th element is the direction of the alternative for the k^th outcome (Accepts "Less" or "Greater" in each entry), (defaults to "Greater" in all entries)
#' @param step the step-size in the subgradient descent (defaults to 100)
#' @param maxIter the maximum number of iterations of subgradient descent (defaults to 1000)
#' @param trueCrit the known, pre-computed chi-bar-squared critical value (defaults to NULL)
#' @param noCorBounds toggles optimistic speed-up using estimated worst-case rho to compute chi-bar-squared quantile if TRUE (defaults to FALSE)
#' @param useNormalQuantile toggles use of critical value from standard Normal distribution, instead of chi-bar-squared if TRUE (defaults to FALSE)
#' @param showDiagnostics toggles diagnostics to be shown
#' @param verbose toggles extra text output
#' @param outputDirName a string for the name of the directory to output the results to (defaults to ""Sensitivity_Analysis_Results")
#' 
#' @return LargestRejectGamma the largest Gamma such that the null hypothesis was rejected
#' 
#' @export

chiBarSquaredTest = function(Q, matchedSetAssignments, treatmentIndicator, 
                             numGamma = 10, alpha = .05, directions = "Greater",
                             step = 100, maxIter = 1000,
                             trueCrit = NULL, noCorBounds = FALSE, useNormalQuantile = FALSE,
                             showDiagnostics = FALSE, verbose = FALSE,
                             outputDirName = "Sensitivity_Analysis_Results")
{
  
  ################################################################################
  #                         Set-up (load packages, format data)
  ################################################################################
  #sources all supporting functions
  allFiles = list.files(path = "./SupportingFunctions")
  for(file in allFiles)
  {
    fileName = paste("./SupportingFunctions/", file, sep = "")
    source(fileName)
  }
  
  #imports all required packages
  require(doParallel)
  require(mvtnorm) #for multivariate normal computations
  require(gurobi) #for using methods from Fogarty & Small 2016
  require(Matrix) #for storing matrices
  require(quadprog) #for solving quadratic programs
  
  K = dim(Q)[1] #the number of variables we observe
  populationSize = dim(Q)[2] #the number of individuals in the population (N in standard notation)
  
  Z = treatmentIndicator #conforming with standard notation
  
  if(all(directions == "Greater")) #If the user enters only "Greater"  as the only direction even when K > 1 the default is that all directions are "Greater"
  {
    directions = rep("Greater", K)
  }
  
  #orders Q, Z, and matchedSetAssignments in increasing stratum number (important for "unlist" functions later on)
  Q = Q[,order(matchedSetAssignments)]
  Z = Z[order(matchedSetAssignments)]
  matchedSetAssignments = matchedSetAssignments[order(matchedSetAssignments)]
  
  #an alternative form of indexing useful for computation.  The t^th element is the list of all individuals in the t^th matched set
  index = list() 
  for(ind in 1:length(unique(matchedSetAssignments)))
  {
    index[[ind]] = which(matchedSetAssignments == ind)
  }
  
  #processing data in case there are multiple treated units and one control in any stratum
  nostratum = length(unique(index)) #the number of strata (matched sets)
  for(i in 1:nostratum)
  {
    ind = index[[i]]
    if(sum(Z[ind]) > 1)
    {
      qsum = base::rowSums(Q[,ind])
      Z[ind] = 1-Z[ind]
      Q[,ind] = qsum - Q[,ind]
      if (showDiagnostics == TRUE)
      {
        cat("Stratum", i, "had multiple treated vs. one control unit.\n")
      }
    }
  }
  
  Q = t(scale(t(Q), center = FALSE, scale = TRUE)) #scales all the rows to have unit standard deviation. does not center them. (helps with numerical stability)
  
  if(verbose == TRUE)
  {
    cat('\nThe data was loaded.\n')
    cat("The population size is:", populationSize, "\n")
    cat("The number of variables is:", K, "\n")
    cat("The hypothesis directions are:\n")
    print(directions)
    cat("The significance level is", alpha, "\n")
  }
  
  ################################################################################
  #               Modify Q to account for different directional choices
  ################################################################################
  directionsScaling = rep(1, length = K)
  directionsScaling[directions == 'Less'] = -1
  print(directionsScaling)
  
  #print(t(head(t(Q))))
  
  Q = diag(directionsScaling)%*%Q #multiplies each row for which there is a "Less" direction by negative 1, then the code proceeds by testing "Greater" directions on this pre-processed data.
  
  
  
  print(t(head(t(Q))))
  
  
  ################################################################################
  #               Compute the univariate test statistics
  ################################################################################
  
  TS = as.vector(Q %*% Z) #the vector of observed test statistics for each of the K outcomes
  if(verbose == TRUE)
  {
    cat("The univariate test statistics after adjusting for directions are:\n")
    print(TS)
  }
  
  ################################################################################
  #               Perform the Double-and-Halve Sensitivity Analysis
  ################################################################################
  InitialGamma = 1
  smallestFailedGamma = Inf
  LargestRejectGamma = 1
  triedGammas = c()
  
  Gamma = InitialGamma
  opt_contrast = rep(NA,K)
  
  while(length(triedGammas) < numGamma)
  {
    if(verbose == TRUE) cat("\n\nTrying Gamma =", Gamma, "\n") #diagnostic
    
    if(Gamma == 1)
    {
      reject = permutationTest(Q = Q, TS = TS, index = index, alpha = alpha, Z=Z, subSampleSize = 500)
    }else{
      reject = computeTestStatistic(Q, TS, index, Gamma, alpha = alpha, Z=Z, step = step,
                                    maxIter = maxIter, trueCrit = trueCrit, noCorBounds = noCorBounds,
                                    useNormalQuantile = useNormalQuantile)
    }
    
    triedGammas = c(triedGammas, Gamma)
    
    if(reject$reject == TRUE) #doubles the Gamma to try next
    {
      LargestRejectGamma = max(Gamma, LargestRejectGamma)
      opt_contrast = reject$contrast
      if(verbose == TRUE) cat("Rejected null at Gamma =", Gamma, "\n")
      if(showDiagnostics == TRUE) cat("The lambdas are: \t", reject$contrast)
    }
    
    if(reject$reject == FALSE) #the next Gamma is the average of the last two Gammas tried (one which rejected and one which failed to reject)
    {
      if(verbose == TRUE) cat("Failed to reject null at Gamma =", Gamma, "\n")
      
      smallestFailedGamma = min(Gamma, smallestFailedGamma)
    }
    
    if (smallestFailedGamma == Inf)
    {
      Gamma = 2*Gamma
    }else if (smallestFailedGamma == 1)
    {
      opt_contrast = reject$contrast
      if(verbose == TRUE) cat("Permutation test at Gamma = 1 failed to reject.\n")
      if(showDiagnostics == TRUE) cat("The lambdas are: \t", reject$contrast)
      break
    }else
    {
      Gamma = (LargestRejectGamma + smallestFailedGamma)/2 #averages the largest tried Gamma which worked with the smallest tried Gamma that failed
    }
    
    #Plots Gamma as a function of the iteration for user interpretation
    plot(x = 1:length(triedGammas),
         y = triedGammas,
         type = 'l',
         xlab = 'Iteration',
         ylab = expression(paste(Gamma)),
         main =  expression(paste(Gamma, " vs. Iteration", sep = '')),
         lwd = 3,
         las = 1
    )
  }
  
  ################################################################################
  #                                 User Output
  ################################################################################
  if(verbose == TRUE)
  {
    if(smallestFailedGamma != 1){
      cat("\n\nSensitivity Analysis Completed:\n\tThe largest Gamma at which the null hypothesis was rejected = ", LargestRejectGamma, "\n\tThe smallest Gamma at which the null failed to be rejected = ", smallestFailedGamma, "\n")
    }else{
      cat("\n\nAt Gamma = 1, the null hypothesis fails to be rejected.\n")
    }
  }
  
  ################################################################################
  #                                 File output
  ################################################################################
  dir.create(outputDirName)
  
  if(smallestFailedGamma != 1){
    write(paste("Sensitivity Analysis Completed:\n\tThe largest Gamma at which the null hypothesis was rejected = ", LargestRejectGamma, "\n\tThe smallest Gamma at which the null failed to be rejected = ", smallestFailedGamma, "\n"), file = paste('./',outputDirName, '/Sensitivity_Analysis_textOutput.pdf', sep = ''))
  }else{
    write("At Gamma = 1, the null hypothesis fails to be rejected.\n", file =  paste('./',outputDirName, '/Sensitivity_Analysis_textOutput.pdf', sep = ''))
  }
  pdf(paste('./',outputDirName, '/Sensitivity_Analysis_imageOutput.pdf', sep = ''))
  plot(x = 1:length(triedGammas),
       y = triedGammas,
       type = 'l',
       xlab = 'Iteration',
       ylab = expression(paste(Gamma)),
       main = expression(paste(Gamma, " vs. Iteration", sep = '')),
       lwd = 3,
       las = 1)
  dev.off()
  
  return(list(LargestRejectGamma=LargestRejectGamma,
              lambda=opt_contrast))
}