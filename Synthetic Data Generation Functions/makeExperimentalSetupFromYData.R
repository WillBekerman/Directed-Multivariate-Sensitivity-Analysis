makeExperimentalSetupFromYData <-function(YData, psi)
{
  s = YData$s_k # Rosenbaums (s_1, s_2) see R16 pg 1459
  
  indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
  
  Z = rep(0, 2 * dim(YData$YData)[1]) 
  Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
  
  Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = length(s))
  
  YOverS = sweep(YData$YData, 2, s, "/")
  
  #print(head(YData))
  #cat(head(apply(YData, MARGIN = 2, FUN = psi)))
  #cat(length(indexRaw))
  
  Q[indexRaw, ] = apply(YOverS, MARGIN = 2, FUN = psi)
  
  # print(head(Q)) # diagnostic
  
  Q[-indexRaw, ] = -Q[indexRaw, ]
  # print(head(Q)) # diagnostic
  
  
  Q = t(Q) # Transposes matrix for formatting's sake
  
  TS = as.vector(Q %*% Z)
  
  
  
  numberOfTreatedIndividuals = length(indexRaw)
  index <- list()
  index[[1]] <- 1:indexRaw[1]
  for(i in 2:length(indexRaw)){
    index[[i]] <- (indexRaw[i-1]+1):indexRaw[i]
  }
  
  return(list(indexRaw = indexRaw,
              index = index,
              Q = Q,
              Z = Z,
              TS = TS))
}