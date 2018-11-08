################################################################################
#             Multiple Comparisons in Sensitivity Analysis
#             This code is directly available at "http://www.mit.edu/~cfogarty/#software"
#               under  "multCompareFunctions.R"
#             It is based upon the methods in the paper:
#                 Fogarty, C. and Small, D. (2016). Sensitivity Analysis for Multiple Comparisons in Matched Observational Studies through Quadratically Constrained Linear Programming. Journal of the American Statistical Association: Theory and Methods,111 (516), 1820-1830 
################################################################################
multipleComparisonsRoot = function(Gamma,index, Q, Z, alpha = 0.05, alternative = "TS", critval = NULL)
{
  
  Z = 1*Z
  
  ns = table(index)
  ms = table(index[Z==1])
  
  nostratum = length(unique(index))
  
  
  
  if(is.null(dim(Q)))
  {
    Q = t(t(Q))
  }
  K = ncol(Q)
  if(any(ms!=1 & ns-ms!=1))
  {
    stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
  }
  if(any(Z!=0 & Z!=1))
  {
    stop("Treatment Vector (Z) Must be Binary")
  }
  for(i in 1:nostratum)
  {
    if(ms[i] > 1)
    {
      ind = which(index==i)
      Z[ind] = 1-Z[ind]
      for(k in 1:K)
      {
        qsum = sum(Q[ind,k])
        
        Q[ind,k] = qsum - Q[ind,k]
      }
    }
  }
  
  treatment = (Z==1)
  PObefore = Q
  K = ncol(PObefore)
  
  sort.new = function(x)
  {
    temp = sort(unique(x))
    new = 1:length(temp)
    ret = rep(0,length(x))
    for(i in new)
    {
      ret[x == temp[i]] = i
    }
    ret
  }
  
  if(length(alternative)==1)
  {
    alternative = rep(alternative, K)
  }
  
  if(any(alternative != "two.sided" & alternative != "greater" & alternative != "less" & alternative != "TS" & alternative != "G"& alternative != "L"))
  {
    stop("Alternative options are two.sided/TS, greater/G, or less/L")  
  }
  alternative[alternative=="less"] = "L"
  alternative[alternative=="greater"] = "G"
  alternative[alternative=="two.sided"] = "TS"
  
  indTS = which(alternative == "TS")
  indG = which(alternative == "G")
  indL = which(alternative == "L")	
  nTS = length(indTS)
  nG = length(indG)
  nL = length(indL)
  orderq = c(indTS, indG, indL)
  PO= PObefore[,orderq, drop = F]
  
  
  Gamma.vec = Gamma
  Reject = rep(0, length(Gamma.vec))
  
  
  
  sds = apply(PO, 2, sd)
  SDS = matrix(sds, nrow(PO), ncol(PO), byrow = T)
  Qmat= (PO/SDS)
  
  ns = table(index)
  ns.types = (ns-1)
  N.total = nrow(Qmat)
  NS = matrix(ns[index], N.total, K)
  nullexpec = colSums(Qmat/NS)
  
  
  treatment = (1*treatment==1)
  Tobs = apply(Qmat[treatment,,drop = F], 2, sum) 
  Tobsmat = matrix(Tobs, N.total, K, byrow = T)  
  
  #bigM = 5*max((Tobs-nullexpec)^2)
  Q = Qmat
  Q2 = Qmat^2
  if(is.null(critval))
  {
    kappa = c(rep(qchisq(1-alpha/K,1),nTS), rep(qchisq(1-2*alpha/K,1),nG+nL))
  }else{kappa = rep(critval, K)}
  kappaMat =  matrix(kappa, nrow(PO), ncol(PO), byrow = T)
  
  Plin = -2*Tobsmat*Q - kappaMat*Q2	
  Plinoneside = -kappaMat*Q2
  U = matrix(0, N.total, length(Gamma))
  RHO = U
  row.ind = rep(0, N.total + (K)*(2*N.total+nostratum+1)+ 4*N.total)
  col.ind = row.ind
  values = row.ind 
  b = rep(0, K*(nostratum+1)+nostratum+2*N.total)
  nvariables = K*(nostratum+1)+N.total+nostratum+(K-nTS)+1
  for(i in 1:nostratum)
  {
    ind = which(index==i)
    row.ind[ind] = rep(i, length(ind))
    col.ind[ind] = ind
    values[ind] = 1
    b[i] = 1
  }
  
  quadcon = vector("list", K)
  
  for(k in 1:K)
  {
    if(k <= nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)  	
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(-Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum+1))
      rowq = N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1))
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims = c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plin[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = -Tobs[k]^2
      
      
      
    }
    
    if(k>nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)    
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum), Tobs[k])
      rowq = c(N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1)-1),N.total + K*(nostratum+1)+1+(k-nTS)) 
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims =  c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plinoneside[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = 0
      
    }
    
    
    
  }
  
  mm = N.total+ K*(2*N.total+nostratum+1)+1
  rr = K*(nostratum+1)+nostratum+1
  cc = N.total + K*(nostratum+1)+2 
  bind = rr
  if(nG!=0)
  {
    genmax  = vector("list", nG)
    
    for(k in (nTS+1):(nTS+nG))
    {
      ktemp = k-nTS
      genmax[[ktemp]]$resvar = cc
      genmax[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmax[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  if(nL!=0)
  {
    genmin = vector("list", nL)
    for(k in (nTS+nG+1):(nTS+nG+nL))
    {
      ktemp = k-nTS-nG
      genmin[[ktemp]]$resvar = cc
      genmin[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmin[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  
  
  mmnext = mm
  ccnext = cc
  rrnext = rr
  
  
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    mm = mmnext
    cc = ccnext
    rr = rrnext
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      for(j in ind)
      {
        row.ind[c(mm, mm+1)] = rep(rr, 2)
        col.ind[c(mm, mm+1)] = c(j, cc)
        values[c(mm, mm+1)] = 	c(1, -Gamma.sens)
        row.ind[c(  mm+2, mm+3)] = rep(rr+1, 2)
        col.ind[c(mm+2, mm+3)] = c(j, cc)
        values[c(mm+2, mm+3)]= c(-1, 1)
        rr = rr+2
        mm = mm+4	
      }
      cc = cc+1	
    }
    
    
    const.dir = c(rep("=", length(b)-2*N.total), rep("<=", 2*N.total))
    model = list()    
    model$A = sparseMatrix(row.ind, col.ind, x=values)
    model$sense = const.dir
    model$quadcon = quadcon
    model$rhs = b	
    model$lb = c(rep(0, N.total), rep(-Inf, K*(1+nostratum)), -1e-4, rep(-Inf, nL+nG), rep(0, nostratum))
    model$obj = c(rep(0, length(model$lb) - nostratum-1-nG-nL), 1, rep(0, nostratum+nG+nL))
    model$ub = c(rep(Inf, N.total), rep(Inf, K*(1+nostratum) + 1), rep(Inf, nL+nG), rep(Inf, nostratum))  
    model$vtypes = c(rep("C", N.total+K*(nostratum+1)+1), rep("C", nL+nG), rep("C", nostratum))
    model$modelsense = "min"
    model$objcons = 0
    
    if(nG>0)
    {
      model$genconmax = genmax
    }
    
    if(nL > 0)
    {
      model$genconmin = genmin
    }
    solm = gurobi(model, params = list(OutputFlag = 0, Cutoff = 1e-4))
    #cat("\t\t\t\t\t\tIn multCompareFunctions.R: The solver found", solm$status, "\n") #diagnostic
    if (solm$status == "INF_OR_UNBD")
    {
      Reject[ee] = -1e-4 # all feasible solutions lead us to fail to reject the null hypothesis
      
    }else if(solm$status == "CUTOFF"){
      Reject[ee] = 1 + 1/(1+Gamma)
      
    }else{
      Reject[ee] = solm$objval
    }
  }
  
  #cat("in multCompareFunctions.R:", Reject, "\n") #diagnostic
  Reject
}







#########
#Same function, just returns the Gamma vector along with whether or not the rejection has been made rather than the objective value
###########
multipleComparisons = function(index, Q, Z, alpha = 0.05, alternative = "TS", Gamma=1)
{
  
  Z = 1*Z
  
  ns = table(index)
  ms = table(index[Z==1])
  
  nostratum = length(unique(index))
  
  
  
  if(is.null(dim(Q)))
  {
    Q = t(t(Q))
  }
  K = ncol(Q)
  if(any(ms!=1 & ns-ms!=1))
  {
    stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
  }
  if(any(Z!=0 & Z!=1))
  {
    stop("Treatment Vector (Z) Must be Binary")
  }
  for(i in 1:nostratum)
  {
    if(ms[i] > 1)
    {
      ind = which(index==i)
      Z[ind] = 1-Z[ind]
      for(k in 1:K)
      {
        qsum = sum(Q[ind,k])
        
        Q[ind,k] = qsum - Q[ind,k]
      }
    }
  }
  
  treatment = (Z==1)
  PObefore = Q
  K = ncol(PObefore)
  
  sort.new = function(x)
  {
    temp = sort(unique(x))
    new = 1:length(temp)
    ret = rep(0,length(x))
    for(i in new)
    {
      ret[x == temp[i]] = i
    }
    ret
  }
  
  if(length(alternative)==1)
  {
    alternative = rep(alternative, K)
  }
  
  if(any(alternative != "two.sided" & alternative != "greater" & alternative != "less" & alternative != "TS" & alternative != "G"& alternative != "L"))
  {
    stop("Alternative options are two.sided/TS, greater/G, or less/L")  
  }
  alternative[alternative=="less"] = "L"
  alternative[alternative=="greater"] = "G"
  alternative[alternative=="two.sided"] = "TS"
  
  indTS = which(alternative == "TS")
  indG = which(alternative == "G")
  indL = which(alternative == "L")	
  nTS = length(indTS)
  nG = length(indG)
  nL = length(indL)
  orderq = c(indTS, indG, indL)
  PO= PObefore[,orderq, drop = F]
  
  
  Gamma.vec = Gamma
  Reject = rep(0, length(Gamma.vec))
  
  
  
  sds = apply(PO, 2, sd)
  SDS = matrix(sds, nrow(PO), ncol(PO), byrow = T)*sqrt(nrow(PO))
  Qmat= (PO/SDS)
  
  ns = table(index)
  ns.types = (ns-1)
  N.total = nrow(Qmat)
  NS = matrix(ns[index], N.total, K)
  nullexpec = colSums(Qmat/NS)
  
  
  treatment = (1*treatment==1)
  Tobs = apply(Qmat[treatment,,drop = F], 2, sum) 
  Tobsmat = matrix(Tobs, N.total, K, byrow = T)  
  
  Q = Qmat
  Q2 = Qmat^2
  
  kappa = c(rep(qchisq(1-alpha/K,1),nTS), rep(qchisq(1-2*alpha/K,1),nG+nL))
  kappaMat =  matrix(kappa, nrow(PO), ncol(PO), byrow = T)
  
  Plin = -2*Tobsmat*Q - kappaMat*Q2	
  Plinoneside = -kappaMat*Q2
  U = matrix(0, N.total, length(Gamma))
  RHO = U
  row.ind = rep(0, N.total + (K)*(2*N.total+nostratum+1)+ 4*N.total)
  col.ind = row.ind
  values = row.ind 
  b = rep(0, K*(nostratum+1)+nostratum+2*N.total)
  nvariables = K*(nostratum+1)+N.total+nostratum+(K-nTS)+1
  for(i in 1:nostratum)
  {
    ind = which(index==i)
    row.ind[ind] = rep(i, length(ind))
    col.ind[ind] = ind
    values[ind] = 1
    b[i] = 1
  }
  
  quadcon = vector("list", K)
  
  for(k in 1:K)
  {
    if(k <= nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)  	
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(-Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum+1))
      rowq = N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1))
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims = c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plin[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = -Tobs[k]^2
      
      
      
    }
    
    if(k>nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)    
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum), Tobs[k])
      rowq = c(N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1)-1),N.total + K*(nostratum+1)+1+(k-nTS)) 
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims =  c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plinoneside[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = 0
      
    }
    
    
    
  }
  
  mm = N.total+ K*(2*N.total+nostratum+1)+1
  rr = K*(nostratum+1)+nostratum+1
  cc = N.total + K*(nostratum+1)+2 
  bind = rr
  if(nG!=0)
  {
    genmax  = vector("list", nG)
    
    for(k in (nTS+1):(nTS+nG))
    {
      ktemp = k-nTS
      genmax[[ktemp]]$resvar = cc
      genmax[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmax[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  if(nL!=0)
  {
    genmin = vector("list", nL)
    for(k in (nTS+nG+1):(nTS+nG+nL))
    {
      ktemp = k-nTS-nG
      genmin[[ktemp]]$resvar = cc
      genmin[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmin[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  
  
  mmnext = mm
  ccnext = cc
  rrnext = rr
  
  
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    mm = mmnext
    cc = ccnext
    rr = rrnext
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      for(j in ind)
      {
        row.ind[c(mm, mm+1)] = rep(rr, 2)
        col.ind[c(mm, mm+1)] = c(j, cc)
        values[c(mm, mm+1)] = 	c(1, -Gamma.sens)
        row.ind[c(  mm+2, mm+3)] = rep(rr+1, 2)
        col.ind[c(mm+2, mm+3)] = c(j, cc)
        values[c(mm+2, mm+3)]= c(-1, 1)
        rr = rr+2
        mm = mm+4	
      }
      cc = cc+1	
    }
    
    
    const.dir = c(rep("=", length(b)-2*N.total), rep("<=", 2*N.total))
    model = list()    
    model$A = sparseMatrix(row.ind, col.ind, x=values)
    model$sense = const.dir
    model$quadcon = quadcon
    model$rhs = b	
    model$lb = c(rep(0, N.total), rep(-Inf, K*(1+nostratum)), -1e-8, rep(-Inf, nL+nG), rep(0, nostratum))
    model$obj = c(rep(0, length(model$lb) - nostratum-1-nG-nL), 1, rep(0, nostratum+nG+nL))
    model$ub = c(rep(Inf, N.total), rep(Inf, K*(1+nostratum) + 1), rep(Inf, nL+nG), rep(Inf, nostratum))  
    model$vtypes = c(rep("C", N.total+K*(nostratum+1)+1), rep("C", nL+nG), rep("C", nostratum))
    model$modelsense = "min"
    
    if(nG>0)
    {
      model$genconmax = genmax
    }
    
    if(nL > 0)
    {
      model$genconmin = genmin
    }
    solm = gurobi(model, params = list(OutputFlag = 0))
    varrho = solm$x[1:N.total]
    RHO[,ee] = varrho
    sums = solm$x[(length(solm$x)-nostratum+1):(length(solm$x))]
    sumsstretch = sums[index]
    uvec = log(varrho/sumsstretch)/log(Gamma.sens)
    U[,ee] = round(uvec, 4)
    Reject[ee] = (solm$objval > 0)
  }
  
  return(list(Gamma=Gamma, Reject = Reject, rho = RHO, u = U))
}


#############
#A function for getting worst case pvalues from the sensitivty analysis for a single outcome (used in the simulation study for the Holm-Bonferroni method)
################ 

sensitivity = function(index, q, Z, alpha = .05, alternative = "two.sided", Gamma.vec=1, calculate.pval = T, continuous.relax = F)
{
  PVAL = calculate.pval
  sdq = sd(q)
  q = q/sd(q)
  Z = 1*Z
  ns = table(index)
  ms = table(index[Z==1])
  ns.types = (ns-1)
  N.total = length(q)
  pval = 0
  nostratum = length(unique(index))
  
  for(i in 1:nostratum)
  {
    ind = which(index==i)
    if(ms[i] > 1)
    {
      qsum = sum(q[ind])
      Z[ind] = 1-Z[ind]
      q[ind] = qsum - q[ind]
    }
  }
  
  
  
  treatment = Z
  
  N.vars = sum((ns-1))
  index = index
  null.expec = sum(q/ns[index])
  Tobs = sum(treatment*q)
  max.e = (sum(treatment*q) > null.expec)
  index.symm = rep(1:nostratum,ns.types)
  
  PM = rep(0, N.vars)
  PV = rep(0, N.vars)
  row.ind = rep(0, 2*N.vars + 1)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nostratum+1)
  for(kk in 1:nostratum)
  {
    row.ind[which(index.symm==kk)] = rep(kk, (ns.types[kk]))
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, ns.types[kk])
    b[kk] = 1
  }
  row.ind[(N.vars+1):(2*N.vars+1)] = rep(nostratum + 1, N.vars+1)
  col.ind[(N.vars+1):(2*N.vars+1)] = 1:(N.vars+1)
  
  opt.expec = rep(0, length(Gamma.vec))
  opt.var = opt.expec
  zscore = opt.expec
  pvalvec = zscore
  Rejectvec = zscore
  Deviate = zscore
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    
    for(kk in 1:nostratum)
    {
      ind = which(index==kk)
      i=kk
      Q = q[ind]
      
      
      
      qi = Q*max.e - Q*(!max.e)
      ord = order(qi)
      qi.sort = sort(qi)
      
      
      mu = rep(0, length(ind)-1)
      sigma2 = rep(0, length(ind)-1)
      
      
      for(j in 1:(length(ind)-1))
      {
        mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
        sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
      }
      mu[abs(mu) < 1e-8] = 0
      sigma2[sigma2 < 1e-8] = 0
      PM[index.symm == kk] = mu*(max.e) - mu*(!max.e)
      PV[index.symm == kk] = (sigma2)
      
      
    }
    
    
    values[(N.vars+1):(2*N.vars+1)] = c(-PM, sign(Tobs))
    b[nostratum+1] = 0
    alpha.opt = alpha
    if(alternative != "two.sided")
    {
      alpha.opt = 2*alpha
    }
    
    
    const.dir = c(rep("=", nostratum+1))
    model = list()
    if(Gamma.sens==1)
    {
      V.test = sum(tapply(PV, index.symm, mean))  
      
      tstat = ((Tobs- null.expec)/sqrt(V.test))
      zed = tstat
      pval = 0
      tstat = zed
      if(alternative == "two.sided")
      {
        pval = 2*pnorm(-abs(tstat))
      }
      if(alternative == "greater")
      {
        pval = 1 - pnorm((tstat))
      }
      if(alternative == "less")
      {
        pval = pnorm((tstat))
      }
      Reject = (pval < alpha)
      
    }
    if(Gamma.sens != 1)
    {
      diff = 10
      kappa = qchisq(1-alpha.opt, 1)
      count=0
      while(diff > 1e-8)
      {
        Plin = -2*Tobs*PM - kappa*PV 
        rowind.q =  1:(N.vars+1)
        colind.q = 1:(N.vars+1)
        values.q = c(rep(0, N.vars),2)
        Q = sparseMatrix(rowind.q, colind.q, x=values.q)
        model$A = sparseMatrix(row.ind, col.ind, x=values)
        model$obj = c(Plin,0)
        model$Q = Q
        model$sense = const.dir
        model$rhs = b
        model$vtype = c(rep("I", N.vars), "C")
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
        model$lb = c(rep(0, N.vars), -Inf)
        
        
        model$modelsense = "min"
        
        
        solm = gurobi(model, params = list(OutputFlag = 0))
        x = solm$x[1:N.vars]
        kappa.new = (Tobs - sum(PM*x))^2/sum(PV*x)
        diff = abs(kappa.new - kappa)
        pval = 0
        if(PVAL == F)
        {
          diff = 0
        }
        kappa = kappa.new
        
      }
      zed = sqrt((Tobs - sum(PM*x))^2/sum(PV*x))
      if(alternative == "less")
      {
        zed = -zed
      }
      zscore[ee] = zed
      tstat = zed
      
      
      if(alternative == "two.sided")
      {
        pval = 2*pnorm(-abs(tstat))
      }
      if(alternative == "greater")
      {
        pval = 1 - pnorm((tstat))
      }
      if(alternative == "less")
      {
        pval = pnorm((tstat))
      }
      Reject = (pval < alpha)
      
      
      
      if(sign(Tobs- sum(PM*x))!=sign(Tobs - null.expec))
      {
        Reject = F
        pval = 0.5
        Deviate = 0
        
        if(alternative == "two.sided")
        {
          pval = 1
        }
      }
      
      if(alternative == "greater" & sum(PM*x) < null.expec)
      {
        pval = .5
      }
      if(alternative == "less" & sum(PM*x) > null.expec)
      {
        pval = .5
      }
      
    }
    pvalvec[ee] = pval
    Rejectvec[ee] = Reject
    Deviate[ee] = abs(zed)
  }    
  
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, Tobs = Tobs*sdq, null.expec = null.expec*sdq, pval = NULL, Deviate = NULL))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec = Gamma.vec, pval = pvalvec, Tobs = Tobs*sdq, null.expec = null.expec*sdq, Deviate = Deviate))
  }
  
}


############
#Aligned Rank
#########
alignedranktest=function(outcome,matchedset,treatment,alternative="two.sided"){
  # Remove units that are not matched
  outcome=outcome[matchedset>0];
  treatment=treatment[matchedset>0];
  matchedset=matchedset[matchedset>0];
  # Compute means in each matched set
  matchedset.mean=tapply(outcome,matchedset,mean);
  # Compute residuals
  matchedset.mean.expand=matchedset.mean[matchedset];
  resids=outcome-matchedset.mean.expand;
  # Rank the residuals
  rankresids=rank(resids);
  # Test statistics = Sum of residuals in treatment group
  teststat=sum(rankresids[treatment==1]);
  # Expected value and variance of test statistic
  mean.matchedset.rankresids=tapply(rankresids,matchedset,mean);
  notreated.matchedset=tapply(treatment,matchedset,sum);
  nocontrol.matchedset=tapply(1-treatment,matchedset,sum);
  no.matchedset=notreated.matchedset+nocontrol.matchedset;
  ev.teststat=sum(mean.matchedset.rankresids*notreated.matchedset);
  mean.matchedset.rankresids.expand=mean.matchedset.rankresids[matchedset];
  rankresids.resid.squared=(rankresids-mean.matchedset.rankresids.expand)^2; 
  squared.resids.sum=tapply(rankresids.resid.squared,matchedset,sum);
  var.teststat=sum(((notreated.matchedset*nocontrol.matchedset)/(no.matchedset*(no.matchedset-1)))*squared.resids.sum);
  
  if(alternative=="two.sided"){
    pval=2*pnorm(-abs((teststat-ev.teststat)/sqrt(var.teststat)));
  }
  if(alternative=="greater"){
    pval=1-pnorm((teststat-ev.teststat)/sqrt(var.teststat));
  }
  if(alternative=="less"){
    pval=pnorm((teststat-ev.teststat)/sqrt(var.teststat));
  }
  pval;
}




###############################
#Define rank based mahalanobis distance
##################################
smahal=
  function(z,X, weight = rep(1, ncol(X))){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k){ X[,j]<-rank(X[,j])}
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    W = diag(1/sqrt(weight))
    icov<-ginv(W%*%cv%*%W)
    
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }
exact.match=function(dmat,z,exact){
  penalty = max(dmat)*100
  adif=abs(outer(exact[z==1],exact[z==0],"-"))
  for(i in 1:nrow(dmat))
  {
    for(j in 1:ncol(dmat))
    {
      if(adif[i,j]!= 0)
      {
        dmat[i,j] = penalty
      }
    }
  }
  dmat
}


########################
#Add a caliper
#######################
addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

#################################
#Standardized Differences for Full Matching
#################################
standardized.diff.func=function(x,treatment,stratum.myindex,missing=rep(0,length(x))){
  xtreated=x[treatment==1 & missing==0];
  xcontrol=x[treatment==0 & missing==0];
  var.xtreated=var(xtreated);
  var.xcontrol=var(xcontrol);
  combinedsd=sqrt(.5*(var.xtreated+var.xcontrol));
  std.diff.before.matching=(mean(xtreated)-mean(xcontrol))/combinedsd;
  nostratum=length(unique(stratum.myindex))-1*min(stratum.myindex==0);
  diff.in.stratum=rep(0,nostratum);
  treated.in.stratum=rep(0,nostratum);
  stratum.size=rep(0,nostratum);
  for(i in 1:nostratum){
    diff.in.stratum[i]=mean(x[stratum.myindex==i & treatment==1 & missing==0])-mean(x[stratum.myindex==i & treatment==0 & missing==0]);
    stratum.size[i] = sum(stratum.myindex==i & missing == 0)
    treated.in.stratum[i]=sum(stratum.myindex==i & treatment==1 & missing==0);
    if(sum(stratum.myindex==i & treatment==0 & missing==0)==0 || sum(stratum.myindex==i & treatment==1 & missing==0)==0){
      treated.in.stratum[i]=0;
      diff.in.stratum[i]=0;
    }
  }
  std.diff.after.matching=(sum(treated.in.stratum*diff.in.stratum)/sum(treated.in.stratum))/combinedsd;
  list(std.diff.before.matching=std.diff.before.matching,std.diff.after.matching=std.diff.after.matching);
}


#################################
#Plot Standardized Differences
############################

plotBalancesign <- function(stdDiff.Before,stdDiff.After,covName,covGroup,titleOfPlot, maxValue = NULL) {
  #in case someone doesn't take the absolute value ahead of time:	
  # Load necessary package
  library(lattice)
  
  # Length of covariate vector
  p = length(covName)
  
  if(is.null(maxValue))
  {
    maxValue = max(stdDiff.Before,stdDiff.After)
  }
  minValue = min(stdDiff.Before, stdDiff.After)
  maxValue = max(abs(maxValue), abs(minValue))
  if(missing(covGroup)) {
    # Setup data frame for lattice package to work
    # In short, it creates a p-by-3 matrix where the rows
    # represent covariates, the first column represents stdDiff before matching
    # the second column represents stdDiff after matching,
    # and the third column represents the names of the covariates. 
    plot.dataframe = data.frame(stdDiff.Before = stdDiff.Before,
                                stdDiff.After = stdDiff.After,
                                covName = covName)
    plot.dataframe$covName = as.factor(plot.dataframe$covName)	
    # This reorder step is necessary to achieve an alphabet-ordering of the covariates on the dot plot
    plot.dataframe$covName = reorder(plot.dataframe$covName,-1*(1:length(plot.dataframe$covName)))
  } else {
    # Set up data frame for lattice package to work
    # We have to reorganize the stdDiff values so that they are grouped 
    # by the covGroup values.
    uniqueCovGroupNames = unique(covGroup)
    nGroups = length(uniqueCovGroupNames)
    
    # We have to add dummy variables for covGroup names
    # Hence, we have p + nGruops, instead of just p
    stdDiff.Before.grouped = rep(0,p + nGroups) 
    stdDiff.After.grouped = rep(0,p + nGroups)
    covName.grouped = rep("",p + nGroups)
    
    index = 1# counter for the for loop
    # Iterate through all the covGroup names
    for(i in 1:length(uniqueCovGroupNames)) {
      # Find covNames that belong to one particular covGroup name
      covInGroup.i = which(covGroup == uniqueCovGroupNames[i])
      
      # Append the name of the covGroup into the stdDiff vectors
      stdDiff.Before.grouped[index] = -3
      stdDiff.After.grouped[index] = -3
      covName.grouped[index] = paste(uniqueCovGroupNames[i],":   ",sep="")
      
      # Copy all the covNames in this particular covGroup into the stdDiff group vector
      stdDiff.Before.grouped[(index+1):(index + 1 + length(covInGroup.i) -1)] = stdDiff.Before[covInGroup.i]
      stdDiff.After.grouped[(index+1):(index + 1 + length(covInGroup.i) - 1)] = stdDiff.After[covInGroup.i]
      covName.grouped[(index+1):(index + 1 + length(covInGroup.i) - 1)] = covName[covInGroup.i]
      
      # Update the index vector
      index = (index + 1 + length(covInGroup.i))
    }
    nPergroup = table(covGroup)
    ln = length(nPergroup)
    index = cbind(c(2, cumsum(nPergroup[1:(ln-1)]) + 1 + (2:ln)),c(cumsum(nPergroup) + (1:ln)))
    
    # Creating the data frame for lattice package. 
    plot.dataframe = data.frame(stdDiff.Before = stdDiff.Before.grouped,
                                stdDiff.After = stdDiff.After.grouped,
                                covName = covName.grouped)
    plot.dataframe$covName = as.factor(plot.dataframe$covName)
    # This reorder step is necessary to achieve an alpha-ordering of the covariates on the dot plot.
    plot.dataframe$covName = reorder(plot.dataframe$covName,-1*(1:length(plot.dataframe$covName)))
    
  }
  dotplot(covName ~ stdDiff.Before,data=plot.dataframe,
          xlim = c(-0.025 - maxValue ,maxValue + 0.025),
          xlab = list("Standardized Differences",cex=0.75),
          main = titleOfPlot,
          panel = function(x,y,...) {
            panel.abline(h = as.numeric(y),lty=2,col="gray") #draws the gray horizonal lines
            
            panel.segments(.2,y[1], .2, y[length(covName)], lty = 6, lwd = 2, col = "red")
            panel.segments(-.2,y[1], -.2, y[length(covName)], lty = 6, lwd = 2, col = "red")
            
            panel.xyplot(x,y,pch=16,col="black",cex=0.75) #plots the before matching stdDiff values
            panel.xyplot(plot.dataframe$stdDiff.After,y,pch=5,col="black",cex=0.6)}, #plots the after matching stdDiff values
          key = list(text = list(c("Before","After"),cex=0.75), #legend
                     points = list(pch = c(16,5),col="black",cex=0.75), #note that pch controls the type of dots
                     space = "right", border=TRUE),
          scales = list(y = list(cex = 0.6)) )
  
} 



########################################
#UnivariateRoot
########################################
univariateRoot = function(Gamma,index, Q, Z, kappa=qnorm(1-alpha/2)^2, alternative = "G")
{
  
  Z = 1*Z
  
  ns = table(index)
  ms = table(index[Z==1])
  
  nostratum = length(unique(index))
  
  
  
  if(is.null(dim(Q)))
  {
    Q = t(t(Q))
  }
  K = ncol(Q)
  if(any(ms!=1 & ns-ms!=1))
  {
    stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
  }
  if(any(Z!=0 & Z!=1))
  {
    stop("Treatment Vector (Z) Must be Binary")
  }
  for(i in 1:nostratum)
  {
    if(ms[i] > 1)
    {
      ind = which(index==i)
      Z[ind] = 1-Z[ind]
      for(k in 1:K)
      {
        qsum = sum(Q[ind,k])
        
        Q[ind,k,drop=F] = qsum - Q[ind,k,drop = F]
      }
    }
  }
  
  treatment = (Z==1)
  PObefore = Q
  K = ncol(PObefore)
  
  sort.new = function(x)
  {
    temp = sort(unique(x))
    new = 1:length(temp)
    ret = rep(0,length(x))
    for(i in new)
    {
      ret[x == temp[i]] = i
    }
    ret
  }
  
  if(length(alternative)==1)
  {
    alternative = rep(alternative, K)
  }
  
  if(any(alternative != "two.sided" & alternative != "greater" & alternative != "less" & alternative != "TS" & alternative != "G"& alternative != "L"))
  {
    stop("Alternative options are two.sided/TS, greater/G, or less/L")  
  }
  alternative[alternative=="less"] = "L"
  alternative[alternative=="greater"] = "G"
  alternative[alternative=="two.sided"] = "TS"
  
  indTS = which(alternative == "TS")
  indG = which(alternative == "G")
  indL = which(alternative == "L")	
  nTS = length(indTS)
  nG = length(indG)
  nL = length(indL)
  orderq = c(indTS, indG, indL)
  PO= PObefore[,orderq, drop = F]
  
  
  Gamma.vec = Gamma
  Reject = rep(0, length(Gamma.vec))
  
  
  
  sds = apply(PO, 2, sd)
  SDS = matrix(sds, nrow(PO), ncol(PO), byrow = T)
  Qmat= (PO/SDS)
  
  ns = table(index)
  ns.types = (ns-1)
  N.total = nrow(Qmat)
  NS = matrix(ns[index], N.total, K)
  nullexpec = colSums(Qmat/NS)
  
  
  treatment = (1*treatment==1)
  Tobs = apply(Qmat[treatment,,drop = F], 2, sum) 
  Tobsmat = matrix(Tobs, N.total, K, byrow = T)  
  
  Q = Qmat
  Q2 = Qmat^2
  
  kappaMat =  matrix(kappa, nrow(PO), ncol(PO), byrow = T)
  
  Plin = -2*Tobsmat*Q - kappaMat*Q2	
  Plinoneside = -kappaMat*Q2
  U = matrix(0, N.total, length(Gamma))
  RHO = U
  row.ind = rep(0, N.total + (K)*(2*N.total+nostratum+1)+ 4*N.total)
  col.ind = row.ind
  values = row.ind 
  b = rep(0, K*(nostratum+1)+nostratum+2*N.total)
  nvariables = K*(nostratum+1)+N.total+nostratum+(K-nTS)+1
  for(i in 1:nostratum)
  {
    ind = which(index==i)
    row.ind[ind] = rep(i, length(ind))
    col.ind[ind] = ind
    values[ind] = 1
    b[i] = 1
  }
  
  quadcon = vector("list", K)
  
  for(k in 1:K)
  {
    if(k <= nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)  	
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(-Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum+1))
      rowq = N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1))
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims = c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plin[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = -Tobs[k]^2
      
      
      
    }
    
    if(k>nTS)
    {
      for(i in 1:nostratum)
      {
        ind = which(index==i)    
        row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
        col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
        values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
      }
      
      row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
      col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
      values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(Q[,k], 1)
      b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum), Tobs[k])
      rowq = c(N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1)-1),N.total + K*(nostratum+1)+1+(k-nTS)) 
      colq = rowq
      valq = rep(1, length(rowq))
      quadcon[[k]] = list()
      quadcon[[k]]$Qc = sparseMatrix(rowq, colq, x = valq, dims =  c(nvariables, nvariables))
      qq= rep(0, nvariables)
      qq[1:N.total] = Plinoneside[,k]
      qq[N.total+K*(nostratum+1)+1] = -1
      
      quadcon[[k]]$q = qq
      quadcon[[k]]$rhs = 0
      
    }
    
    
    
  }
  
  mm = N.total+ K*(2*N.total+nostratum+1)+1
  rr = K*(nostratum+1)+nostratum+1
  cc = N.total + K*(nostratum+1)+2 
  bind = rr
  if(nG!=0)
  {
    genmax  = vector("list", nG)
    
    for(k in (nTS+1):(nTS+nG))
    {
      ktemp = k-nTS
      genmax[[ktemp]]$resvar = cc
      genmax[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmax[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  if(nL!=0)
  {
    genmin = vector("list", nL)
    for(k in (nTS+nG+1):(nTS+nG+nL))
    {
      ktemp = k-nTS-nG
      genmin[[ktemp]]$resvar = cc
      genmin[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
      genmin[[ktemp]]$con = 0
      cc=cc+1
    }
  }
  
  
  mmnext = mm
  ccnext = cc
  rrnext = rr
  
  
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    mm = mmnext
    cc = ccnext
    rr = rrnext
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      for(j in ind)
      {
        row.ind[c(mm, mm+1)] = rep(rr, 2)
        col.ind[c(mm, mm+1)] = c(j, cc)
        values[c(mm, mm+1)] = 	c(1, -Gamma.sens)
        row.ind[c(  mm+2, mm+3)] = rep(rr+1, 2)
        col.ind[c(mm+2, mm+3)] = c(j, cc)
        values[c(mm+2, mm+3)]= c(-1, 1)
        rr = rr+2
        mm = mm+4	
      }
      cc = cc+1	
    }
    
    
    const.dir = c(rep("=", length(b)-2*N.total), rep("<=", 2*N.total))
    model = list()    
    model$A = sparseMatrix(row.ind, col.ind, x=values)
    model$sense = const.dir
    model$quadcon = quadcon
    model$rhs = b	
    model$lb = c(rep(0, N.total), rep(-Inf, K*(1+nostratum)), -1e-4, rep(-Inf, nL+nG), rep(0, nostratum))
    model$obj = c(rep(0, length(model$lb) - nostratum-1-nG-nL), 1, rep(0, nostratum+nG+nL))
    model$ub = c(rep(Inf, N.total), rep(Inf, K*(1+nostratum) + 1), rep(Inf, nL+nG), rep(Inf, nostratum))  
    model$vtypes = c(rep("C", N.total+K*(nostratum+1)+1), rep("C", nL+nG), rep("C", nostratum))
    model$modelsense = "min"
    model$objcons = 0
    
    if(nG>0)
    {
      model$genconmax = genmax
    }
    
    if(nL > 0)
    {
      model$genconmin = genmin
    }
    solm = gurobi(model, params = list(OutputFlag = 0, Cutoff = 1e-4))
    #cat("\t\t\t\t\t\tIn multCompareFunctions.R: The solver found", solm$status, "\n") #diagnostic
    if (solm$status == "INF_OR_UNBD")
    {
      Reject[ee] = -1e-4 # all feasible solutions lead us to fail to reject the null hypothesis
      
    }else if(solm$status == "CUTOFF"){
      Reject[ee] = 1 + 1/(1+Gamma)
      
    }else{
      Reject[ee] = solm$objval
    }
  }
  
  #cat("in multCompareFunctions.R:", Reject, "\n") #diagnostic
  Reject
}

########################################
#initialPoint
########################################
initialPoint = function(index, Q, Z, alpha = 0.05, alternative = "TS", Gamma=1)
{
  Z = 1*Z
  
  ns = table(index)
  ms = table(index[Z==1])
  
  nostratum = length(unique(index))
  
  
  
  if(is.null(dim(Q)))
  {
    Q = t(t(Q))
  }
  K = ncol(Q)
  if(any(ms!=1 & ns-ms!=1))
  {
    stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
  }
  if(any(Z!=0 & Z!=1))
  {
    stop("Treatment Vector (Z) Must be Binary")
  }
  for(i in 1:nostratum)
  {
    if(ms[i] > 1)
    {
      ind = which(index==i)
      Z[ind] = 1-Z[ind]
      for(k in 1:K)
      {
        qsum = sum(Q[ind,k])
        
        Q[ind,k,drop=F] = qsum - Q[ind,k,drop = F]
      }
    }
  }
  
  treatment = (Z==1)
  PObefore = Q
  K = ncol(PObefore)
  
  sort.new = function(x)
  {
    temp = sort(unique(x))
    new = 1:length(temp)
    ret = rep(0,length(x))
    for(i in new)
    {
      ret[x == temp[i]] = i
    }
    ret
  }
  
  if(length(alternative)==1)
  {
    alternative = rep(alternative, K)
  }
  
  if(any(alternative != "two.sided" & alternative != "greater" & alternative != "less" & alternative != "TS" & alternative != "G"& alternative != "L"))
  {
    stop("Alternative options are two.sided/TS, greater/G, or less/L")  
  }
  
  PO= PObefore[,,drop = F]
  
  
  Gamma.vec = Gamma
  Reject = rep(0, length(Gamma.vec))
  
  
  
  sds = apply(PO, 2, sd)
  SDS = matrix(sds, nrow(PO), ncol(PO), byrow = T)
  Qmat= (PO/SDS)
  
  ns = table(index)
  ns.types = (ns-1)
  N.total = nrow(Qmat)
  NS = matrix(ns[index], N.total, K)
  nullexpec = colSums(Qmat/NS)
  
  
  treatment = (1*treatment==1)
  
  
  Plin = Q	
  U = matrix(0, N.total, length(Gamma))
  RHO = U
  S = matrix(0, nostratum, length(Gamma))
  row.ind = rep(0, N.total + 2*nostratum)
  col.ind = row.ind
  values = row.ind 
  b = rep(0, nostratum+2*N.total)
  for(i in 1:nostratum)
  {
    ind = which(index==i)
    row.ind[ind] = rep(i, length(ind))
    col.ind[ind] = ind
    values[ind] = 1
    b[i] = 1
  }
  
  
  
  mmnext = N.total+1
  ccnext = N.total+1
  rrnext = nostratum+1
  
  
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    mm = mmnext
    cc = ccnext
    rr = rrnext
    for(i in 1:nostratum)
    {
      ind = which(index == i)
      for(j in ind)
      {
        row.ind[c(mm, mm+1)] = rep(rr, 2)
        col.ind[c(mm, mm+1)] = c(j, cc)
        values[c(mm, mm+1)] = 	c(1, -Gamma.sens)
        row.ind[c(  mm+2, mm+3)] = rep(rr+1, 2)
        col.ind[c(mm+2, mm+3)] = c(j, cc)
        values[c(mm+2, mm+3)]= c(-1, 1)
        rr = rr+2
        mm = mm+4	
      }
      cc = cc+1	
    }
    
    
    const.dir = c(rep("=", nostratum), rep("<=", 2*N.total))
    model = list()    
    model$A = sparseMatrix(row.ind, col.ind, x=values)
    model$sense = const.dir
    model$rhs = b	
    model$lb = c(rep(0, N.total), rep(0, nostratum))
    model$obj = c(Q, rep(0, nostratum))
    model$ub = c(rep(Inf, N.total+nostratum))  
    model$vtypes = c(rep("C", N.total+nostratum))
    model$modelsense = "max"
    model$objcons = 0
    solm = gurobi(model, params = list(OutputFlag = 0))
    varrho = solm$x[1:N.total]
    RHO[,ee] = varrho
    sums = solm$x[(length(solm$x)-nostratum+1):(length(solm$x))]
    S[,ee] = sums
    sumsstretch = sums[index]
    uvec = log(varrho/sumsstretch)/log(Gamma.sens)
    U[,ee] = round(uvec, 4)
    Reject[ee] = (solm$objval > 0)
  }
  
  return(list(Gamma=Gamma, Reject = Reject, rho = as.vector(RHO), s = as.vector(S), u = U))
}


