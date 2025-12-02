################################################################################
#  Get functions from Cohen et al. (2020) paper,
#  only modifications are for speed-ups which are explicitly noted
################################################################################
library("httr")
library("jsonlite")

source_gitfile <- function(filename){
  url = paste0("https://raw.githubusercontent.com/WillBekerman/Directed-",
               "Multivariate-Sensitivity-Analysis/master/",filename)
  script = GET(url = url)
  script = content(script,"text")
  eval(parse(text = script), envir = .GlobalEnv)
}
resp <- GET(
  paste0("https://api.github.com/repos/WillBekerman/Directed-Multivariate-",
  "Sensitivity-Analysis/contents/SupportingFunctions?ref=master"),
  add_headers(
    "Accept" = "application/vnd.github.v3+json",
    "User-Agent" = "httr"
  )
)
items <- fromJSON(content(resp, as="text", encoding="UTF-8"))
# source chiBarSquaredTest.R and all .R files in ./SupportingFunctions
file_vec <- c("chiBarSquaredTest.R", "Synthetic%20Data%20Generation%20Functions/makeExperimentalSetupFromYData.R", items$path)
getfile <- sapply(file_vec, source_gitfile)




################################################################################
#  Get Table (1−alpha)-quantiles of some relevant distributions at alpha = 0.05
################################################################################
library("mvtnorm")
library("fourPNO")
library("quadprog")

normal_quantiles = chisq_quantiles = chibarsq_zero_corr_quantiles = 
  chibarsq_pos_corr_quantiles = numeric()
alpha = 0.05
pos_corr = 0.2
for (p in c(2,5,25,100,150)){
  normal_quantiles <- c(normal_quantiles, qnorm(1-alpha))
  chisq_quantiles <- c(chisq_quantiles, sqrt(qchisq(1-alpha, p)))
  ### chibarsq quantile: zero correlation
  wts = numeric(p+1L)
  wts[1] = pmvnorm(rep(0,p), rep(Inf,p), sigma=diag(p))[[1]]
  wts[p+1L] = pmvnorm(rep(0,p), rep(Inf,p), sigma=diag(p))[[1]]
  for(i in seq(1L, p-1L, by=1L)) wts[i+1] = (2^(-p)*factorial(p))/
    (factorial(i)*factorial(p-i))
  chibarsq_zero_corr_quantiles = c(chibarsq_zero_corr_quantiles,
                                   sqrt(qchibarsq(1-alpha, diag(p), wts)))
  ## chibarsq quantile: correlation given by pos_corr
  wts = numeric(p+1L)
  V = (1 - pos_corr) * diag(p) + pos_corr * matrix(1, p, p)
  if (p <= 12) { # `wchibarsq()` starts becoming intractable around p=12
    wts = wchibarsq(solve(V))
  } else { # simulate weights
    nsim = 1000000
    for(sim in 1:nsim){
      y = rmvnorm(n=1, mu=rep(0,p), sigma = V)
      matprod = solve.QP( Dmat = solve(V),
                          dvec = y%*%solve(V),
                          Amat = diag(1,p),
                          bvec = rep(0,p) )$solution
      numpos = sum(matprod>1E-6)
      wts[numpos+1] <- wts[numpos+1] + 1
    }
    wtsvec=wts/nsim
    wts<-rev(wtsvec) #w_i(p,V) = w_{p-i}(p,V^{-1})
  }
  chibarsq_pos_corr_quantiles <- c(chibarsq_pos_corr_quantiles,
                                   sqrt(qchibarsq(1-alpha, solve(V), wts)))
}
normal_quantiles
chisq_quantiles
chibarsq_zero_corr_quantiles
chibarsq_pos_corr_quantiles




################################################################################
#  Define our main functions
################################################################################
# Main function on planning stage -- get optimal contrast at sensitivity value
main_planning <- function(dat, directions, method='PGD', psi.f='Huber', trim.in.vec=NULL, numGamma=15, alpha=0.05){
  
  # check dat, directions, method, numGamma, alpha arguments
  stopifnot(!is.null(dat))
  stopifnot(!is.null(directions) && (length(directions)==ncol(dat) && sum(directions %in% c('Greater','Less'))==length(directions)))
  stopifnot(!is.null(method) && (method=='L-BFGS-B')||(method=='PGD'))
  stopifnot(!is.null(numGamma) && numGamma>=1)
  stopifnot(!is.null(alpha) && alpha>0 && alpha<1)
  
  # define M-statistic
  stopifnot(is.null(psi.f) || psi.f %in% c('Huber','InnerTrim'))
  stopifnot(is.null(trim.in.vec) || (length(trim.in.vec)==2 && trim.in.vec[1]>=trim.in.vec[2] && trim.in.vec[2]>=0))
  if (is.null(psi.f)){ # we must have input trim and inner
    trim <- trim.in.vec[1]
    inner <- trim.in.vec[2]
  } else { # we must have input psi.f \in {'Huber','InnerTrim'}
    if (psi.f == 'Huber'){
      trim = 2.5; inner = 0
    } else if (psi.f == 'InnerTrim'){
      trim = 2.5; inner = 0.5
    }
  }

  # define internal functions
  get_Q <- function(dat){
    psi <- Vectorize( FUN = function(y){
      return(sign(y)*(trim / (trim - inner)) * max(c(0, min(c(abs(y), trim)) - inner)))
    })
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    # ii = 0
    # while (any(s==0)) {
    #   ii = ii+1
    #   s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) quantile(x,0.5+ii/15,na.rm=T))
    # }
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    Q[indexRaw, ] = apply(YOverS, MARGIN = 2, FUN = psi)
    Q[-indexRaw, ] = -Q[indexRaw, ]
    Q = t(Q) # Transposes matrix for formatting's sake
    numberOfTreatedIndividuals = length(indexRaw)
    index <- list()
    for(i in 1:length(indexRaw)){
      index[[i]] <- c(indexRaw[i]-1,indexRaw[i])
    }
    matchedSetAssignments = rep(0, ncol(Q))
    for(ind in 1:length(index))
    {
      matchedSetAssignments[unlist(index[ind])] = ind
    }
    list(Q=Q, matchedSetAssignments=matchedSetAssignments, Z=Z)
  }
  get_parts_of_deviate <- function(Q,gamma,Lweights){
    tvec <- apply(Q,1,function(x)sum(x[seq(2,ncol(Q),2)])) # since we set first individual in each pair to be the control indiv.
    tstat <- sum( Lweights*apply(Q,1,function(x)sum(x[seq(2,ncol(Q),2)])) )
    newQ <- Lweights%*%Q
    ymat=matrix(newQ,nrow=nostratum,ncol=2,byrow=T)
    rhomat <- t(apply(ymat, 1, function(x) sort(x,index.return=T)$ix))
    m=2;j=1;
    pr <- c(rep(1, j), rep(gamma, m - j))/(j + ((m - j)*gamma))
    rhomat[rhomat==1] <- pr[1]
    rhomat[rhomat==2] <- pr[2]
    
    ### each variable
    tvec <- evec <- rep(NA,nrow(Q))
    covmat <- matrix(NA,nrow(Q),nrow(Q))
    for (varnum in 1:nrow(Q)){
      ymat=matrix(Q[varnum,],nrow=nostratum,ncol=2,byrow=T)
      mu=diag(ymat %*% t(rhomat))
      tstat <- as.vector(sum(ymat[, 2])) # since we set first individual in each pair to be the control indiv.
      expect <- sum(mu)
      tvec[varnum] <- tstat
      evec[varnum] <- expect
      ### each pair of variables
      for (varnumnew in 1:varnum){
        ymat.varnum=ymat
        ymat.varnumnew=matrix(Q[varnumnew,],nrow=nostratum,ncol=2,byrow=T)
        sigma2 <- as.vector(diag((ymat.varnum * ymat.varnumnew) %*% t(rhomat))) - (mu * as.vector(diag(ymat.varnumnew %*% t(rhomat))))
        vartotal <- sum(sigma2)
        covmat[varnumnew,varnum] <- vartotal
      }
    }
    covmat[lower.tri(covmat)] <- t(covmat)[lower.tri(covmat)]
    list(t=tvec,e=evec,covmat=covmat)
  }
  comparisonUnew <- function(Q,gamma,lam){
    
    partsofdev <- get_parts_of_deviate(Q,gamma,lam)
    t <- partsofdev$t
    e <- partsofdev$e
    covmat <- partsofdev$covmat
    dev <- as.numeric( t(lam)%*%(t-e) / (sqrt( t(lam)%*%covmat%*%lam )) )
    
    list(deviate=dev)
  }
  optim_comparison <- function(Q,gamma,lambda){
    tmp<-comparisonUnew(Q,gamma,lambda)
    return(tmp$deviate)
  }

  # get Q matrix and scale by specified directions
  dat <- as.matrix(dat)
  nvar <- ncol(dat)
  nostratum <- nrow(dat)
  dat <- as.numeric(dat) # change to vector formatting
  Qres <- get_Q(dat)
  Q <- Qres$Q
  matchedSetAssignments <- Qres$matchedSetAssignments
  Z <- Qres$Z
  Q = t(scale(t(Q), center = FALSE, scale = TRUE)) # scale rows to have unit SD, don't center them (help with numerical stability)
  directionsScaling = rep(1,length=nvar)
  directionsScaling[directions=='Less'] = -1
  Q = sweep(Q, MARGIN = 1, STATS = directionsScaling, FUN = "*")
  
  # get optimal contrast and planning sample sensitivity value
  triedGammas = c()
  opt_contrast = rep(NA,nvar)
  gamma = 1
  LargestRejectGamma = 1
  smallestFailedGamma = Inf
  p = Inf
  while(length(triedGammas) < numGamma){
    triedGammas = c(triedGammas, gamma)
    # get lambda that maximizes deviate
    if (method == 'L-BFGS-B'){ # box constraints, Byrd et al. (1995)
      opt_results<-optim(rnorm(nvar,.5,.5),optim_comparison,Q=Q,gamma=gamma,
                         method = method,
                         lower=c( rep(0,nvar) ), upper=c( rep(1,nvar) ),
                         control=list(fnscale=-1))
      contrast<-opt_results$par/sum(opt_results$par)
      contrast<-contrast*directionsScaling
    }
    else if (method == 'PGD'){ # projected sub-gradient method of Cohen et al. (2020)
      index = list()
      for(ind in 1:length(unique(matchedSetAssignments))){
        index[[ind]] = which(matchedSetAssignments == ind)
      }
      opt_results<-computeTestStatistic(Q, as.vector(Q%*%Z), index, gamma, Z, alpha, step=5,
                                        maxIter=1000, noCorBounds=F, useNormalQuantile=T)
      contrast<-opt_results$contrast/sum(opt_results$contrast)
      contrast<-contrast*directionsScaling
    }
    sa_result <- comparisonUnew(Q,gamma,lam=contrast/directionsScaling)
    dev<-sa_result$deviate
    p <- 1-pnorm(dev) # use Gaussian reference dist., comment on this in paper
    rej <- p <= alpha
    if(rej){ # doubles the gamma to try next
      LargestRejectGamma = max(gamma, LargestRejectGamma)
      opt_contrast <- contrast
      devreturn <- dev
    } else{ # next gamma is avg of last two Gammas tried (one rejected, one FTR)
      smallestFailedGamma = min(gamma, smallestFailedGamma)
    }
    if (smallestFailedGamma == Inf){
      gamma = 2*gamma
    } else if (smallestFailedGamma == 1){ # gamma = 1 failed to reject
      opt_contrast <- contrast
      devreturn <- dev
      break
    } else{
      gamma = (LargestRejectGamma + smallestFailedGamma)/2 # averages the largest tried gamma which worked with the smallest tried gamma that failed
    }
  }
  lambda <- opt_contrast
  sensval <- LargestRejectGamma
  return(list(lambda=lambda,
              sensval=sensval))
}

# Main function on analysis stage -- get sensitivity value or pvalue at Gamma
main_fixedlam <- function(dat, lam, psi.f='Huber', trim.in.vec=NULL, Gamma=NULL, numGamma=NULL, alpha=NULL){
  
  # check dat, lam, Gamma, numGamma, alpha arguments
  stopifnot(!is.null(dat))
  stopifnot(!is.null(lam) && length(lam)==ncol(dat))
  stopifnot((!is.null(Gamma) && Gamma>=1) || (!is.null(numGamma) && numGamma>=1 && !is.null(alpha) && alpha>0 && alpha<1))
  
  # define M-statistic
  stopifnot(is.null(psi.f) || psi.f %in% c('Huber','InnerTrim'))
  stopifnot(is.null(trim.in.vec) || (length(trim.in.vec)==2 && trim.in.vec[1]>=trim.in.vec[2] && trim.in.vec[2]>=0))
  if (is.null(psi.f)){ # we must have input trim and inner
    trim <- trim.in.vec[1]
    inner <- trim.in.vec[2]
  } else { # we must have input psi.f \in {'Huber','InnerTrim'}
    if (psi.f == 'Huber'){
      trim = 2.5; inner = 0
    } else if (psi.f == 'InnerTrim'){
      trim = 2.5; inner = 0.5
    }
  }
  
  # define internal function
  get_Q <- function(dat){
    psi <- Vectorize( FUN = function(y){
      return(sign(y)*(trim / (trim - inner)) * max(c(0, min(c(abs(y), trim)) - inner)))
    })
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    # ii = 0
    # while (any(s==0)) {
    #   ii = ii+1
    #   s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) quantile(x,0.5+ii/15,na.rm=T))
    # }
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    Q[indexRaw, ] = apply(YOverS, MARGIN = 2, FUN = psi)
    Q[-indexRaw, ] = -Q[indexRaw, ]
    Q = t(Q) # Transposes matrix for formatting's sake
    numberOfTreatedIndividuals = length(indexRaw)
    index <- list()
    for(i in 1:length(indexRaw)){
      index[[i]] <- c(indexRaw[i]-1,indexRaw[i])
    }
    matchedSetAssignments = rep(0, ncol(Q))
    for(ind in 1:length(index))
    {
      matchedSetAssignments[unlist(index[ind])] = ind
    }
    list(Q=Q, matchedSetAssignments=matchedSetAssignments)
  }

  # get Q matrix, multiply by lam
  dat <- as.matrix(dat)
  nvar <- ncol(dat)
  nostratum <- nrow(dat)
  dat <- as.numeric(dat) # change to vector formatting
  Q <- get_Q(dat)$Q
  Q = t(scale(t(Q), center = FALSE, scale = TRUE)) # scale rows to have unit SD, don't center them (help with numerical stability)
  Q = lam %*% Q

  # get pval at Gamma or sensitivity value on analysis sample
  if (!is.null(Gamma)){
    # get p-value by modifying sensitivitymv::mscorev()
    #-1/2*matrix(Q,ncol=2,byrow=T) for pairs
    ms <- Q
    ms <- matrix(ms, ncol=2, byrow = TRUE)
    ms <- ms[,2:1,drop = FALSE] # mscorev wants (T,C) not (C,T)
    n <- dim(ms)[1]
    m <- dim(ms)[2]
    ms <- array(ms, c(n, m, m - 1))
    ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
    ni <- apply(!is.na(ms), 1, sum)
    use <- (ni >= 2) & (!is.na(ms[, 1]))
    ms <- ms[use, ]
    ni <- ni[use]
    ms <- ms/outer(ni, rep(1, m), "*")
    res = sensitivitymv::separable1v(ms, gamma = Gamma)
    dev <- res$deviate
    p <- 1-pnorm(dev)
    return(list(pval=p))
  } else {
    triedGammas = c()
    gamma = 1
    LargestRejectGamma = 1
    smallestFailedGamma = Inf
    p = Inf
    while(length(triedGammas) < numGamma){
      triedGammas = c(triedGammas, gamma)
      # get sensitivity value by modifying sensitivitymv::mscorev()
      ms <- Q
      ms <- matrix(ms, ncol=2, byrow = TRUE)
      ms <- ms[,2:1,drop = FALSE] # mscorev wants (T,C) not (C,T)
      n <- dim(ms)[1]
      m <- dim(ms)[2]
      ms <- array(ms, c(n, m, m - 1))
      ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
      ni <- apply(!is.na(ms), 1, sum)
      use <- (ni >= 2) & (!is.na(ms[, 1]))
      ms <- ms[use, ]
      ni <- ni[use]
      ms <- ms/outer(ni, rep(1, m), "*")
      res = sensitivitymv::separable1v(ms, gamma = gamma)
      dev <- res$deviate
      p <- 1-pnorm(dev)
      rej <- p <= alpha
      if(rej){ # doubles the gamma to try next
        LargestRejectGamma = max(gamma, LargestRejectGamma)
        devreturn <- dev
      } else{ # next gamma is avg of last two Gammas tried (one rejected, one FTR)
        smallestFailedGamma = min(gamma, smallestFailedGamma)
      }
      if (smallestFailedGamma == Inf){
        gamma = 2*gamma
      } else if (smallestFailedGamma == 1){ # gamma = 1 failed to reject
        devreturn <- dev
        break
      } else{
        gamma = (LargestRejectGamma + smallestFailedGamma)/2 # averages the largest tried gamma which worked with the smallest tried gamma that failed
      }
    }
    sensval <- LargestRejectGamma
    return(list(sensval=sensval))
  }
}




################################################################################
#  Run Simulations
################################################################################
library("tidyverse")

run_sim <- function(nostratum,
                    Taus,
                    directions,
                    correlation, # the known, constant correlation between outcome variables
                    psi.f = 'Huber',
                    noise.dist = 'Normal',
                    run_comparators = TRUE, # Whole w/ Search, Whole w/ DS
                    planning_sample_prop = 0.25,
                    nsim = 1000,
                    seed = 0,
                    numGamma = 15,
                    alpha = 0.05,
                    parallel = TRUE){
  
  # set seed for reproducibility
  set.seed(seed)
  
  # simulation metrics to return
  wholewithsearch_sv <- wholewithds_sv <- split_sv <- numeric(length=nsim)

  # define M-statistic
  stopifnot(psi.f %in% c('Huber','InnerTrim'))
  if (psi.f == 'Huber'){
    trim = 2.5; inner = 0
  }
  if (psi.f == 'InnerTrim'){
    trim = 2.5; inner = 0.5
  }
  
  # generate data with parametric noise distribution
  stopifnot(noise.dist %in% c('Normal'))
  generateData <- function(rho, tauvec, nostratum){
    errors = fourPNO::rmvnorm(n=nostratum, mu=rep(0,length(tauvec)),
                              sigma=rho+diag(x=1-rho,
                                             nrow=length(tauvec),
                                             ncol=length(tauvec)))
    ymat = matrix(NA, nrow=nostratum, ncol=length(tauvec))
    for (j in 1:length(tauvec)){
      ymat[,j] <- errors[,j] + tauvec[j]
    }
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:length(tauvec),paste0) )
    normalData = list(YData = data.frame(Ypairs = Ypairs_normal))
    return(list(normalData = normalData))
  }
  
  # get design sensitivity for given DGP
  psi <- Vectorize( FUN = function(y){
    return(sign(y)*(trim / (trim - inner)) * max(c(0, min(c(abs(y), trim)) - inner)))
  })
  syntheticData_asymp = generateData(rho = correlation,
                                     tauvec = Taus,
                                     nostratum = 1e7)
  syntheticData_asymp_small = generateData(rho = correlation,
                                           tauvec = Taus,
                                           nostratum = 1e5)
  omega_k <- apply(abs(syntheticData_asymp$normalData$YData),
                   MARGIN = 2, FUN = median)
  YOverOmega = sweep(syntheticData_asymp_small$normalData$YData, MARGIN = 2,
                     omega_k, "/")
  psi_YOverOmega = apply(YOverOmega, MARGIN = 2, FUN = psi)
  psi_YOverOmega[,directions == 'Less']=-psi_YOverOmega[,directions == 'Less']
  eta <- function(lam,psi_res){
    mean(sqrt( (lam%*%t(psi_res))^2 ))
  }
  theta <- function(lam,psi_res){
    mean( lam%*%t(psi_res) )
  }
  get_design_sens <- function(lam,psi_res){
    (eta(lam,psi_res)+theta(lam,psi_res)) /
      (eta(lam,psi_res)-theta(lam,psi_res))
  }
  lam_optdesignsens <- try({
    optim(par = rep(0.01, ncol(psi_YOverOmega)), fn = get_design_sens,
          psi_res = psi_YOverOmega, method = "L-BFGS-B", lower = 0, upper = 1,
          control = list(lmm = 5, pgtol = 1e-8, fnscale = -1e-3))$par
  })
  # run if first optimization doesn't work, should help convergence
  if (class(lam_optdesignsens) == "try-error") {
    lam_optdesignsens <-
      optim(par = rep(0.01, ncol(psi_YOverOmega)), fn = get_design_sens,
            psi_res = psi_YOverOmega, method = "L-BFGS-B", lower = 0, upper = 1,
            control = list(lmm = 5, pgtol = 1e-8, fnscale = -1))$par
  }
  lam_optdesignsens <- lam_optdesignsens/sum(lam_optdesignsens)
  # scale by pre-specified directions
  directionsScaling = rep(1,length=length(Taus))
  directionsScaling[directions=='Less'] = -1
  lam_optdesignsens = lam_optdesignsens*directionsScaling
  
  # function to get chibarsq critical value if we run Cohen et al. (2020) as comparator
  get_chibarsqcrit <- function(trueCor, K){ # can either use `wchibarsq()` for small K, simulate wts, or closed form for UB
    wts<-numeric(K+1L)
    if (trueCor == 0){ # trueCor is zero (i.e., covariance is identity); closed-form soln
      V=diag(K)
      wts[1]<-pmvnorm(rep(0,K),rep(Inf,K),sigma=solve(V))[[1]]
      wts[K+1L]<-pmvnorm(rep(0,K), rep(Inf,K),sigma=V)[[1]]
      for(i in seq(1L, K-1L, by=1L)) wts[i+1]= (2^(-K)*factorial(K))/(factorial(i)*factorial(K-i))
    }
    else{
      V = (1 - trueCor) * diag(K) + trueCor * matrix(1, K, K)
      if (K<=12) { # `wchibarsq()` starts becoming intractable around K=12
        wts = wchibarsq(solve(V))
      } else { # have to simulate weights
        nsim=1000000 # hyperparameter, typically works well in practice
        for(sim in 1:nsim){
          y = fourPNO::rmvnorm(n=1, mu=rep(0,K), sigma = V)
          matprod=quadprog::solve.QP( Dmat = solve(V),
                                      dvec = y%*%solve(V),
                                      Amat = diag(1,K),
                                      bvec = rep(0,K) )$solution
          numpos=sum(matprod>1E-6)
          wts[numpos+1] <- wts[numpos+1] + 1
        }
        wtsvec=wts/nsim
        wts<-rev(wtsvec) #w_i(p,V) = w_{p-i}(p,V^{-1})
      }
    }
    crit = qchibarsq(1-alpha, solve(V), wts) # note, if wts argument non-null, V vs. solve(V) doesnt matter
    return(crit)
  }
  if (run_comparators){
    chibarsqcrit <- get_chibarsqcrit(correlation, length(Taus))
  }
  
  #-----loop over simulations-----
  if (parallel){
    library(foreach)
    library(doParallel)
    library(doRNG)
    num_cores <- detectCores() - 1 # used 10 cores for K=35 to avoid computer crashing, also used 10 cores for I = 10k
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    # registerDoRNG(seed = seed)
    
    objs <- ls(envir = .GlobalEnv)
    funs <- objs[sapply(objs, function(x) is.function(get(x, envir = .GlobalEnv)))]

    results <- foreach(sim = 1:nsim, .combine = cbind, .packages = c("quadprog"),
                       .export = funs#c("psi","main_planning","main_fixedlam","chiBarSquaredTest",sub(".R.*", "", sub(".*/", "", file_vec)))
                       ) %dopar% {
                         
      showDiagnostics = verbose = FALSE
                         

      # create and split synthetic data into planning and analysis samples
      syntheticData = generateData(rho = correlation,
                                   tauvec = Taus,
                                   nostratum = nostratum)
      Data_whole = syntheticData$normalData
      planning_sample_size_sets <- floor(planning_sample_prop*nostratum)
      ix_planning <- sort(sample(x=1:nostratum,size=planning_sample_size_sets))
      ix_analysis <- (1:nostratum)[-ix_planning]
      # split data and run our sensitivity analysis
      # planning sample
      Data_planning <- Data_whole$YData[ix_planning,]
      planning_result <- main_planning(Data_planning, directions, method='PGD', 
                                       psi.f=psi.f, numGamma=numGamma)
      # analysis sample
      Data_analysis <- Data_whole$YData[ix_analysis,]
      analysis_result <- main_fixedlam(Data_analysis,
                                       planning_result$lambda,
                                       psi.f=psi.f,
                                       numGamma=numGamma,
                                       alpha=alpha)
      split_sv <- analysis_result$sensval
      if (run_comparators){ # run alternative sensitivity analyses
        # Comparator: whole with search (Cohen et al. 2020)
        Data_whole$s_k = apply(abs(Data_whole$YData), MARGIN = 2, FUN = median)
        psi <- Vectorize( FUN = function(y){
          return(sign(y)*(trim / (trim - inner)) * max(c(0, min(c(abs(y), trim)) - inner)))
        })
        experimentalSetup_whole = makeExperimentalSetupFromYData(Data_whole,psi)
        Q_whole = experimentalSetup_whole$Q
        Z_whole = experimentalSetup_whole$Z
        index_whole = experimentalSetup_whole$index
        matchedSetAssignments_whole = rep(0, ncol(Q_whole))
        for(ind in 1:length(index_whole)){
          matchedSetAssignments_whole[unlist(index_whole[ind])] = ind
        }
        wholewithsearch_result = chiBarSquaredTest(Q = Q_whole,
                                                   matchedSetAssignments = matchedSetAssignments_whole,
                                                   treatmentIndicator = Z_whole,
                                                   numGamma = numGamma,
                                                   alpha = alpha,
                                                   directions = directions,
                                                   step = 5,
                                                   maxIter = 1000,
                                                   trueCrit = chibarsqcrit,
                                                   showDiagnostics = showDiagnostics,
                                                   verbose = verbose,
                                                   outputDirName = "Sims_Results_WholewithSearch")
        # Comparator: whole with design sensitivity
        wholewithds_result = main_fixedlam(Data_whole$YData,
                                           lam_optdesignsens,
                                           psi.f=psi.f,
                                           numGamma=numGamma,
                                           alpha=alpha)
        wholewithsearch_sv <- wholewithsearch_result$LargestRejectGamma
        wholewithds_sv <- wholewithds_result$sensval
        return(c(split_sv, wholewithsearch_sv, wholewithds_sv))
      } else{
        return(split_sv)
      }
    }
    
    if (run_comparators){
      split_sv <- as.numeric(results[1,])
      wholewithsearch_sv <-as.numeric(results[2,])
      wholewithds_sv <- as.numeric(results[3,])
    } else {
      split_sv <- as.numeric(results[1,])
    }
    
    stopCluster(cl)
    
  } else{
    for (sim in 1:nsim){
      
      cat('\n\n\n\n\nSIMULATION NUMBER: ', sim, '\n\n\n\n\n\n')
      
      # create and split synthetic data into planning and analysis samples
      syntheticData = generateData(rho = correlation,
                                   tauvec = Taus,
                                   nostratum = nostratum)
      Data_whole = syntheticData$normalData
      planning_sample_size_sets <- floor(planning_sample_prop*nostratum)
      ix_planning <- sort(sample(x=1:nostratum,size=planning_sample_size_sets))
      ix_analysis <- (1:nostratum)[-ix_planning]
      # split data and run our sensitivity analysis
      # planning sample
      Data_planning <- Data_whole$YData[ix_planning,]
      planning_result <- main_planning(Data_planning, directions, method='PGD', 
                                       psi.f=psi.f, numGamma=numGamma)
      # analysis sample
      Data_analysis <- Data_whole$YData[ix_analysis,]
      analysis_result <- main_fixedlam(Data_analysis,
                                       planning_result$lambda,
                                       psi.f=psi.f,
                                       numGamma=numGamma,
                                       alpha=alpha)
      if (run_comparators){ # run alternative sensitivity analyses
        # Comparator: whole with search (Cohen et al. 2020)
        Data_whole$s_k = apply(abs(Data_whole$YData), MARGIN = 2, FUN = median)
        experimentalSetup_whole = makeExperimentalSetupFromYData(Data_whole,psi)
        Q_whole = experimentalSetup_whole$Q
        Z_whole = experimentalSetup_whole$Z
        index_whole = experimentalSetup_whole$index
        matchedSetAssignments_whole = rep(0, ncol(Q_whole))
        for(ind in 1:length(index_whole)){
          matchedSetAssignments_whole[unlist(index_whole[ind])] = ind
        }
        wholewithsearch_result = chiBarSquaredTest(Q = Q_whole,
                                                   matchedSetAssignments = matchedSetAssignments_whole,
                                                   treatmentIndicator = Z_whole,
                                                   numGamma = numGamma,
                                                   alpha = alpha,
                                                   directions = directions,
                                                   step = 5,
                                                   maxIter = 1000,
                                                   trueCrit = chibarsqcrit,
                                                   showDiagnostics = F,
                                                   verbose = F,
                                                   outputDirName = "Sims_Results_WholewithSearch")
        # Comparator: whole with design sensitivity
        wholewithds_result = main_fixedlam(Data_whole$YData,
                                           lam_optdesignsens,
                                           psi.f=psi.f,
                                           numGamma=numGamma,
                                           alpha=alpha)
        wholewithsearch_sv[sim] <- wholewithsearch_result$LargestRejectGamma
        wholewithds_sv[sim] <- wholewithds_result$sensval
      }
      split_sv[sim] <- analysis_result$sensval
    }
  }
  
  
  
  
  # return metrics as a list called res
  res=list()
  res$split_sv <- split_sv
  res$wholewithsearch_sv <- wholewithsearch_sv
  res$wholewithds_sv <- wholewithds_sv
  res=res[unlist(lapply(res,function(x)!all(x==0)))] # remove NA entries
  return(res)
}
#----- Get simulation results -----


# setting mixed: Taus = rep(c(-.25,-.5,.1,.25,.5), 3)
# setting sparse (old): Taus = rep(c(-.1,-.1,.1,.5,.5), 3)

run_comparators = TRUE
if (run_comparators){
  methodnames <- c('Analysis w/ Planning','Whole w/ Search', 'Whole w/ DS')
} else {
  methodnames <- c('Analysis w/ Planning')
}
nsim = 1000
showDiagnostics = verbose = TRUE
notify_done <- function(msg = "✅ R finished running!") {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Windows") {
    system(paste('msg *', shQuote(msg)))
  } else if (sys == "Darwin") {
    system(paste0("osascript -e 'display notification ", shQuote(msg), " with title \"R Alert\"'"))
  } else if (sys == "Linux") {
    system(paste('notify-send', shQuote(msg)))
  } else {
    message(msg)
  }
}



nostratum = 300 # I \in {300,1000}

for(numrep in c(1,3,5)){ # K \in {5,15,25}

  for (corr_val in c(-0.2,0,0.2)){ # rho \in {-0.2, 0, 0.2}
    
    cat('\n\n\n\n\nNEW SIMULATION SETTING: ', '\nK =',numrep*5, 'Outcomes\nCorrelation =',corr_val,'\n\n\n\n\n\n')
    
    Taus = rep(c(.1,.1,.1,.1,.5), numrep)
    directions = rep(c('Greater','Greater','Greater','Greater','Greater'), numrep)
    
    sim_result=run_sim(nostratum=nostratum,
                       Taus=Taus,
                       directions=directions,
                       correlation=corr_val,
                       run_comparators=run_comparators,
                       nsim=nsim)
    
    save(sim_result, file=paste0('sparse-K',numrep*5,'-I',nostratum,'-cor',corr_val,'.RData'))
    notify_done("🎯 All simulations completed successfully!")
  }
  
}


plots_list <- vector(mode="list",length=6)
ix = 1
useBW = T

if (useBW) colorvals <- rep("grey20",3) else colorvals <- c("#8CBAD8", "#D3A46E", "#B8A3CD")

for(corr_val in c(-0.2,0,0.2)){ # rho \in {-0.2, 0, 0.2}
  
  for (numrep in c(1,3,5)){ # K \in {5,15,25}
    
    if (numrep == 1) gammas_vector <- seq(from=1,to=4,by=.25)
    if (numrep == 3) gammas_vector <- seq(from=1.5,to=9.5,by=.25)
    if (numrep == 5) gammas_vector <- seq(from=2,to=16,by=.25)
    
    load(file=paste0('sparse-K',numrep*5,'-I',nostratum,'-cor',corr_val,'.RData'))
    power_result <- lapply(sim_result,function(x)colSums(outer(x,gammas_vector,`>`))/nsim)
    df <- data.frame(Gamma=rep(gammas_vector,length(methodnames)),
                     Sensitivity=unlist(power_result),
                     Method=rep(methodnames,each=length(gammas_vector)))
    plt <- ggplot(df, aes(x=Gamma, y=Sensitivity)) +
      geom_line(aes(linetype=Method, col=Method), size=1.3, alpha = 0.95) +
      scale_color_manual(values = colorvals)+
      scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
      labs(
        x = expression(Gamma),
        y = "Power"
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        aspect.ratio = 1,
        axis.line = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "grey50", linewidth = 0.2),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "grey20"),
        axis.ticks.length = unit(.1, "in")
      ) + 
      scale_x_continuous(breaks = scales::breaks_pretty(4))
    plots_list[[ix]] <- plt
    ix = ix+1
  }
  
}

combined_plt <- cowplot::plot_grid(plotlist = plots_list, ncol = 3)
ggsave(paste0(if(useBW)'bw-','combinedplt-sparse-I',nostratum,'.pdf'), combined_plt, width = 12, height = 9, device = "pdf")


