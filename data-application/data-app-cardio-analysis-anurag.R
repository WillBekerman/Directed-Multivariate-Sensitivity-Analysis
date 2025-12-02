
################################################################################
#  Set parameters and define functions
################################################################################
library(tidyverse)
library(mvtnorm)
library(fourPNO)
library(quadprog)
library(globaltest)
library(sensitivitymv)
library(ggplot2)
library(reshape)
library(httr)
library(jsonlite)

# Set parameters
planning_sample_prop = 0.25
alpha = 0.05; gamma = 1
showDiagnostics = verbose=T # show output when main functions are being run

# Inheritance procedure functions (Goeman & Finos, 2012) using family law-inspired
# heir function, node weights as the number of leaves below or at each node,
# and Shaffer improvement.
# Code directly from Goeman's Github:
inheritance_procedure <- function( ps,
                                   weights,
                                   parents,
                                   children,
                                   offspring,
                                   Shaffer = TRUE,
                                   homogeneous = FALSE) {
  
  # how many tests? 
  m <- length(ps)
  nms <- names(ps)
  weights <- weights[nms]
  
  # find top and leaves
  top <- sapply(parents[nms], length) == 0
  names(top) <- nms
  leaf <- sapply(children[nms], length) == 0
  names(leaf) <- nms
  
  # find parents of leaves
  leaf.parents <- children[sapply(children, function(ch) any(ch %in% nms[leaf]))]
  leaf.siblings <- lapply(leaf.parents, function(lp) {
    nleaves <- sum(leaf[lp])
    if (nleaves > 1) lp else lp[!leaf[lp]]
  })
  
  # initialize
  basealpha <- rep(0,m)
  names(basealpha) <- nms
  basealpha[top] <- weights[top]/sum(weights[top])   # basealpha should add up to 1. Is multiplied by alpha later
  shaffer <- rep(1, m)
  rejected <- rep(FALSE, m)
  names(rejected) <- nms
  extinct <- weights == 0    # nodes with weight zero are never inherited to
  names(extinct) <- nms
  adjp <- rep(1, m)
  names(adjp) <- nms
  
  # a flag to do plain Meinshausen
  Meinshausen <- FALSE
  
  # start the procedure
  alpha <- 0
  ready <- FALSE
  while (!ready) {
    
    # phase 1: reject
    testalpha <- basealpha * shaffer
    newly.rejected <- (ps/testalpha <= alpha) & (!rejected)
    newly.rejected[ps == 0 & testalpha == 0] <- FALSE     # do not reject when 0/0
    
    if (any(newly.rejected)) {
      adjp[newly.rejected] <- alpha
      rejected <- rejected | newly.rejected
      
      # phase 2: recalculate extinctness
      extinct[(!extinct) & rejected] <- sapply((1:m)[(!extinct) & rejected], function(i) 
        leaf[i] || all(rejected[offspring[[nms[i]]]]) || sum(weights[offspring[[nms[[i]]]]]) == 0)
      
      # phase 3: recalculate Shaffer
      if (Shaffer) {
        shaffer <- rep(1, m)
        names(shaffer) <- names(ps)
        if (homogeneous) {
          bottom <- setdiff(nms[rejected], unlist(parents[nms[rejected]]))
          bottom <- setdiff(bottom, nms[leaf])
          minweights <- unlist(lapply(bottom, function(bt) {
            min(weights[intersect(offspring[[bt]], nms[leaf])])
          }))
          sumweights <- sum(weights[leaf & !rejected])
          if (sumweights > 0) 
            shaffer[1:m] <- sumweights / (sumweights - sum(minweights))
          else
            shaffer[1:m] <- 1
          for (bt in bottom) {
            if (length(children[[bt]]) > 1) {
              mw <- unlist(lapply(children[[bt]], function(ch) {
                min(weights[intersect(c(ch,offspring[[ch]]), nms[leaf])])
              }))  
              names(mw) <- children[[bt]]
              smallest <- names(which.min(mw))
              second <- min(mw[names(mw) != smallest])
              shaffer[smallest] <- sumweights / (sumweights - sum(minweights) + mw[smallest] - second)
            } else
              shaffer[children[[bt]]] <- Inf
          }
        } else {
          for (lp in names(leaf.parents)[rejected[names(leaf.parents)]]) {
            desc <- leaf.parents[[lp]]
            lsibs <- leaf.siblings[[lp]]
            if (!any(rejected[desc])) {
              lvs <- desc[leaf[desc]]
              if (length(desc) == 1)
                shaffer[desc] <- Inf
              else {             
                smallest <- names(which.min(weights[lvs])) 
                shaffer[lsibs] <- sum(weights[desc]) /  (sum(weights[desc]) - weights[smallest])
                if (smallest %in% lsibs) {
                  second <- min(weights[setdiff(lvs, smallest)])
                  shaffer[smallest] <- sum(weights[desc]) / (sum(weights[desc]) - second)
                }  
              }
            }
          }
        }
      }
      
      # phase 4: inherit
      while (any(rejected & (basealpha > 0))) {
        for (i in 1:m) {
          if (rejected[i] && basealpha[i] > 0) {
            
            # find heirs
            if (extinct[i])
              if (top[i])         # whole tree is rejected
                heirs <- nms[top & (!extinct)]
            else                # heirs of leaf nodes or nodes under which all leaves are rejected
              if (Meinshausen)
                heirs <- character(0)
              else if (homogeneous) {
                heirs <- nms[(basealpha > 0) & (1:m != i)]
              }  
              else
                heirs <- parents[[nms[i]]]
              else {                 # heirs of internal nodes
                heirs <- children[[nms[i]]]
                heirs <- heirs[!extinct[heirs]]
              }  
              
              # inherit basealpha to heirs
              basealpha[heirs] <- basealpha[heirs] + basealpha[i] * weights[heirs]/sum(weights[heirs])
              basealpha[i] <- 0
          }
        }
      }
    } else {        # no new rejections: increase alpha
      if (all(rejected))
        ready <- TRUE
      else {
        next.rejection <- ps / (basealpha * shaffer)
        next.rejection[ps == 0 & basealpha == 0] <- Inf
        alpha <- min(next.rejection)
        ready <- (alpha >= 1)
        newalpha <- TRUE
      }
    }
  }
  adjp[!rejected] <- 1
  return(adjp)
}
# Wrapper function to simplify input format
run_inheritance <- function(nodes, 
                            wts = NULL,
                            Shaffer = TRUE,
                            homogeneous = FALSE) {
  
  # Validate node structure
  if (!all(c("parents", "children", "p") %in% unique(unlist(lapply(nodes, names))))) {
    stop("All nodes must contain 'parents', 'children', and 'p' components")
  }
  
  # Extract components
  node_names <- names(nodes)
  ps <- sapply(nodes, function(n) n$p)
  names(ps) <- node_names
  
  # Get parents and children
  parents <- lapply(nodes, function(n) n$parents)
  children <- lapply(nodes, function(n) n$children)
  
  # Build offspring list (node + all descendants)
  build_offspring <- function(name) {
    unique(c(name, unlist(lapply(children[[name]], build_offspring))))
  }
  offspring <- setNames(lapply(node_names, build_offspring), node_names)
  
  # Calculate weights if not given
  if (is.null(wts)) {
    is_leaf <- sapply(children, function(x) length(x) == 0)
    wts <- sapply(offspring, function(desc) sum(is_leaf[desc])) # define node weights as the number of leaves below or at each node
    names(wts) <- node_names
  }
  
  # Convert to required formats
  parents <- setNames(lapply(parents, function(p) {
    if(length(p) == 0) character(0) else p
  }), node_names)
  
  children <- setNames(lapply(children, function(ch) {
    if(length(ch) == 0) character(0) else ch
  }), node_names)
  
  # Run inheritance procedure
  inheritance_procedure(
    ps = ps,
    weights = wts,
    parents = parents,
    children = children,
    offspring = offspring,
    Shaffer = Shaffer,
    homogeneous = homogeneous
  )
}



################################################################################
#  Define our main functions
################################################################################
# Load auxillary functions directly from our Github
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

# Note: `function_general()` allows for running our method with composites
# To do so, set dat as a list, with each component containing $data and $lambda
# to make a composite, or just $data.
# Also, we allow for searching over statistics and applying different test stats for each outcome
# Main function on planning stage -- get optimal contrast at sensitivity value
main_planning_general <- function(dat, directions, method='PGD', psi.f='Huber', trim.in.vec=NULL, trim.in.list=NULL, numGamma=15, alpha=0.05){

  # check dat, method, numGamma, alpha arguments
  stopifnot(!is.null(dat))
  stopifnot(!is.null(method) && (method=='L-BFGS-B')||(method=='PGD'))
  stopifnot(!is.null(numGamma) && numGamma>=1)
  stopifnot(!is.null(alpha) && alpha>0 && alpha<1)
  
  # define M-statistic
  if (!is.null(trim.in.list)){ # many statistics
    pars <- param_list <- trim.in.list
  } 
  else { # one statistic
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
    if (is.list(dat)){
      ncomposites = length(dat)
      pars <- lapply(seq_len(ncomposites), function(i) {
        list(c(trim, inner)
        )
      })
    } else{
      param_list <- list(c(trim, inner))
    }
  }
  
  # define internal functions
  get_Q1 <- function(dat){
    psi <- function(y, trim, inner) {
      sign(y) * (trim / (trim - inner)) * pmax(0, pmin(abs(y), trim) - inner)
    }
    
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    tempmat = mapply(
      FUN = function(col, params) {
        psi(col, trim = params[1], inner = params[2])
      },
      col    = as.data.frame(YOverS),  # columns as list
      params = param_list,
      SIMPLIFY = FALSE
    )
    Q[indexRaw, ] = do.call(cbind, tempmat)
    
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
  get_Q2 <- function(dat,pars){
    psi <- function(y, trim, inner) {
      sign(y) * (trim / (trim - inner)) * pmax(0, pmin(abs(y), trim) - inner)
    }
    
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    tempmat = mapply(
      FUN = function(col, params) {
        psi(col, trim = params[1], inner = params[2])
      },
      col    = as.data.frame(YOverS),  # columns as list
      params = pars,
      SIMPLIFY = FALSE
    )
    Q[indexRaw, ] = do.call(cbind, tempmat)
    
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
  get_Q <- function(dat, pars = NULL) {
    if (is.null(pars))
      return(get_Q1(dat))
    else
      return(get_Q2(dat, pars))
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
  
  # get Q matrix and scale by specified directions: general case allowing for composites
  if (is.list(dat)){ # build composite measure(s)
    ncomposites = length(dat)
    Qlist = vector("list", length=ncomposites)
    for (compositenum in 1:ncomposites){
      nvar <- ncol(as.matrix(dat[[compositenum]]$data))
      nostratum <- nrow(as.matrix(dat[[compositenum]]$data))
      Qres <- get_Q(as.numeric(as.matrix(dat[[compositenum]]$data)), pars[compositenum][[1]])
      if (is.null(dat[[compositenum]]$lambda)){ # usual case
        Qlist[[compositenum]] <- Qres$Q
      } else { # will produce univariate composite measure
        Qlist[[compositenum]] <- dat[[compositenum]]$lambda %*% Qres$Q
      }
    }
    Q <- do.call(rbind, Qlist)
  } else{ # usual case
    nvar <- ncol(dat)
    nostratum <- nrow(dat)
    dat <- as.matrix(dat)
    Qres <- get_Q(as.numeric(dat))
    Q <- Qres$Q
  }
  
  nvar <- nrow(Q)
  nostratum <- ncol(Q)/2
  
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
main_fixedlam_general <- function(dat, lam, psi.f='Huber', trim.in.vec=NULL, trim.in.list=NULL, Gamma=NULL, numGamma=NULL, alpha=NULL){
  
  # check dat, Gamma, numGamma, alpha arguments
  stopifnot(!is.null(dat))
  stopifnot((!is.null(Gamma) && Gamma>=1) || (!is.null(numGamma) && numGamma>=1 && !is.null(alpha) && alpha>0 && alpha<1))
  
  # define M-statistic
  if (!is.null(trim.in.list)){ # many statistics
    pars <- param_list <- trim.in.list
  } 
  else { # one statistic
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
    if (is.list(dat)){
      ncomposites = length(dat)
      pars <- lapply(seq_len(ncomposites), function(i) {
        list(c(trim, inner)
        )
      })
    } else{
      param_list <- list(c(trim, inner))
    }
  }
  
  # define internal function
  get_Q1 <- function(dat){
    psi <- function(y, trim, inner) {
      sign(y) * (trim / (trim - inner)) * pmax(0, pmin(abs(y), trim) - inner)
    }
    
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    tempmat = mapply(
      FUN = function(col, params) {
        psi(col, trim = params[1], inner = params[2])
      },
      col    = as.data.frame(YOverS),  # columns as list
      params = param_list,
      SIMPLIFY = FALSE
    )
    Q[indexRaw, ] = do.call(cbind, tempmat)
    
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
  get_Q2 <- function(dat,pars){
    psi <- function(y, trim, inner) {
      sign(y) * (trim / (trim - inner)) * pmax(0, pmin(abs(y), trim) - inner)
    }
    
    ymat=matrix(dat,nrow=nostratum,ncol=nvar)
    Ypairs_normal = data.frame(ymat)
    names(Ypairs_normal) <- as.character( outer('Y',1:nvar,paste0) )
    YData = list(YData = data.frame(Ypairs = Ypairs_normal))
    s = apply(abs(YData$YData), MARGIN = 2, FUN = function(x) median(x,na.rm=T))
    indexRaw = seq(2, 2 * dim(YData$YData)[1], by = 2) # We set the first individual in each pair to be the control indiv.
    Z = rep(0, 2 * dim(YData$YData)[1]) 
    Z[indexRaw] = 1 # Sets up the indicator of treatment in accordance with the treatment indices.
    Q = matrix(0, nrow = 2 * dim(YData$YData)[1], ncol = nvar)
    YOverS = sweep(YData$YData, 2, s, "/")
    tempmat = mapply(
      FUN = function(col, params) {
        psi(col, trim = params[1], inner = params[2])
      },
      col    = as.data.frame(YOverS),  # columns as list
      params = pars,
      SIMPLIFY = FALSE
    )
    Q[indexRaw, ] = do.call(cbind, tempmat)
    
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
  get_Q <- function(dat, pars = NULL) {
    if (is.null(pars))
      return(get_Q1(dat))
    else
      return(get_Q2(dat, pars))
  }
  
  # get Q matrix and multiply by lam: general case allowing for composites
  if (is.list(dat)){ # build composite measure(s)
    ncomposites = length(dat)
    Qlist = vector("list", length=ncomposites)
    for (compositenum in 1:ncomposites){
      nvar <- ncol(as.matrix(dat[[compositenum]]$data))
      nostratum <- nrow(as.matrix(dat[[compositenum]]$data))
      Qres <- get_Q(as.numeric(as.matrix(dat[[compositenum]]$data)), pars[compositenum][[1]])
      if (is.null(dat[[compositenum]]$lambda)){ # usual case
        Qlist[[compositenum]] <- Qres$Q
      } else { # will produce univariate composite measure
        Qlist[[compositenum]] <- dat[[compositenum]]$lambda %*% Qres$Q
      }
    }
    Q <- do.call(rbind, Qlist)
  } else{ # usual case
    nvar <- ncol(dat)
    nostratum <- nrow(dat)
    dat <- as.matrix(dat)
    Qres <- get_Q(as.numeric(dat))
    Q <- Qres$Q
  }
  
  nvar <- nrow(Q)
  nostratum <- ncol(Q)/2
  
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
      if(is.na(rej)){ # edge case where all difs = 0
        LargestRejectGamma=smallestFailedGamma=1
        break
      }
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

# Used to write our results to JSON for plotting
as_path <- function(name) {
  m <- regexec("^Ages([0-9]+to[0-9]+)_([MF])_(.+)$", name)
  hits <- regmatches(name, m)[[1]]
  if (length(hits) == 4) {
    rng <- gsub("to", "-", hits[2])
    sex <- hits[3]
    list("All Outcomes", paste0("Ages ", rng), paste0(sex, ": ", rng), name)
  } else {
    list("All Outcomes", name)  # fallback for non-standard names
  }
}



################################################################################
#  Process data and divide into pilot/analysis samples
################################################################################
diet_outcomes_names <- c('HEI_TotalFruit', 'HEI_WholeFruit', 'HEI_TotalVeg',
                         'HEI_GreensBeans', 'HEI_WholeGrain', 'HEI_Dairy',
                         'HEI_TotalProt', 'HEI_SeaPlantProt', 'HEI_FattyAcid',
                         'HEI_RefinedGrain', 'HEI_Sodium', 'HEI_AddSug',
                         'HEI_SatFat', 'HEI_Total')
all_outcome_names <- c('BMXBMI',
                       'waist_height',
                       'moderate_leisure',
                       'vigorous_leisure',
                       'LBXCOT',
                       'sys_bp',
                       'LBXGH',
                       'non_HDL_chol',
                       'eGFR',
                       diet_outcomes_names)
outcome_names_8to11 <- c('BMXBMI',
                         'waist_height',
                         'LBXCOT',
                         'sys_bp',
                         'non_HDL_chol',
                         diet_outcomes_names)

dat_outcomes_8to11_m <- subset(dat, select=outcome_names_8to11)
treated_8to11_m <- treated_idx[treated_idx %in% which(dat$AGEGRP==1 & dat$RIAGENDR==1)]
control_8to11_m <- control_idx[control_idx %in% which(dat$AGEGRP==1 & dat$RIAGENDR==1)]
dat_outcomes_12to17_m <- subset(dat, select=all_outcome_names)
treated_12to17_m <- treated_idx[treated_idx %in% which(dat$AGEGRP==2 & dat$RIAGENDR==1)]
control_12to17_m <- control_idx[control_idx %in% which(dat$AGEGRP==2 & dat$RIAGENDR==1)]
dat_outcomes_8to11_f <- subset(dat, select=outcome_names_8to11)
treated_8to11_f <- treated_idx[treated_idx %in% which(dat$AGEGRP==1 & dat$RIAGENDR==2)]
control_8to11_f <- control_idx[control_idx %in% which(dat$AGEGRP==1 & dat$RIAGENDR==2)]
dat_outcomes_12to17_f <- subset(dat, select=all_outcome_names)
treated_12to17_f <- treated_idx[treated_idx %in% which(dat$AGEGRP==2 & dat$RIAGENDR==2)]
control_12to17_f <- control_idx[control_idx %in% which(dat$AGEGRP==2 & dat$RIAGENDR==2)]

# split sample into pilot/analysis samples
planning_num_sets_8to11_m <- floor(planning_sample_prop*length(treated_8to11_m))
ix_planning_sets_8to11_m <- sample(x=1:length(treated_8to11_m),size=planning_num_sets_8to11_m)
ix_analysis_sets_8to11_m <- (1:length(treated_8to11_m))[-ix_planning_sets_8to11_m]
planning_num_sets_12to17_m <- floor(planning_sample_prop*length(treated_12to17_m))
ix_planning_sets_12to17_m <- sample(x=1:length(treated_12to17_m),size=planning_num_sets_12to17_m)
ix_analysis_sets_12to17_m <- (1:length(treated_12to17_m))[-ix_planning_sets_12to17_m]
planning_num_sets_8to11_f <- floor(planning_sample_prop*length(treated_8to11_f))
ix_planning_sets_8to11_f <- sample(x=1:length(treated_8to11_f),size=planning_num_sets_8to11_f)
ix_analysis_sets_8to11_f <- (1:length(treated_8to11_f))[-ix_planning_sets_8to11_f]
planning_num_sets_12to17_f <- floor(planning_sample_prop*length(treated_12to17_f))
ix_planning_sets_12to17_f <- sample(x=1:length(treated_12to17_f),size=planning_num_sets_12to17_f)
ix_analysis_sets_12to17_f <- (1:length(treated_12to17_f))[-ix_planning_sets_12to17_f]

# pilot sample paired differences
mpdifs_pilot_8to11_m <- dat_outcomes_8to11_m[treated_8to11_m[ix_planning_sets_8to11_m],]-dat_outcomes_8to11_m[control_8to11_m[ix_planning_sets_8to11_m],]
mpdifs_pilot_8to11_f <- dat_outcomes_8to11_f[treated_8to11_f[ix_planning_sets_8to11_f],]-dat_outcomes_8to11_f[control_8to11_f[ix_planning_sets_8to11_f],]
mpdifs_pilot_12to17_m <- dat_outcomes_12to17_m[treated_12to17_m[ix_planning_sets_12to17_m],]-dat_outcomes_12to17_m[control_12to17_m[ix_planning_sets_12to17_m],]
mpdifs_pilot_12to17_f <- dat_outcomes_12to17_f[treated_12to17_f[ix_planning_sets_12to17_f],]-dat_outcomes_12to17_f[control_12to17_f[ix_planning_sets_12to17_f],]
mpdifs_pilot_8to11 <- rbind(mpdifs_pilot_8to11_m, mpdifs_pilot_8to11_f)
mpdifs_pilot_12to17 <- rbind(mpdifs_pilot_12to17_m, mpdifs_pilot_12to17_f)



################################################################################
#  Exploratory data analysis on pilot sample
################################################################################
# Plots of paired differences: Ages 8-11
plist_8to11_m <- plist_8to11_f <- vector(mode="list",length=length(outcome_names_8to11))
for (colnum in 1:length(outcome_names_8to11)){
  
  cat(paste('Outcome:',outcome_names_8to11[colnum]),'\n')
  
  df_m <- melt(data.frame(MatchedDifs=
                            as.numeric(unlist(dat_outcomes_8to11_m[treated_8to11_m[ix_planning_sets_8to11_m],colnum]))-
                            as.numeric(unlist(dat_outcomes_8to11_m[control_8to11_m[ix_planning_sets_8to11_m],colnum]))))
  
  df_f <- melt(data.frame(MatchedDifs=
                            as.numeric(unlist(dat_outcomes_8to11_f[treated_8to11_f[ix_planning_sets_8to11_f],colnum]))-
                            as.numeric(unlist(dat_outcomes_8to11_f[control_8to11_f[ix_planning_sets_8to11_f],colnum]))))
  
  plist_8to11_m[[colnum]] <- ggplot(df_m, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() + theme_minimal() + theme(legend.position="none") +
    labs(title=paste(outcome_names_8to11[colnum]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)
  
  plist_8to11_f[[colnum]] <- ggplot(df_f, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() + theme_minimal() + theme(legend.position="none") +
    labs(title=paste(outcome_names_8to11[colnum]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)

  # Two-sided p-values
  cat(paste('Male pvalue:',wilcox.test( df_m$value )$p.value,'\n'))
  cat(paste('Female pvalue:',wilcox.test( df_f$value )$p.value,'\n\n\n'))
  
}

# Also plot cotinine levels with fewer outliers
df_m <- melt(data.frame(MatchedDifs=
                          as.numeric(unlist(dat_outcomes_8to11_m[treated_8to11_m[ix_planning_sets_8to11_m],3]))-
                          as.numeric(unlist(dat_outcomes_8to11_m[control_8to11_m[ix_planning_sets_8to11_m],3]))))
df_f <- melt(data.frame(MatchedDifs=
                          as.numeric(unlist(dat_outcomes_8to11_f[treated_8to11_f[ix_planning_sets_8to11_f],3]))-
                          as.numeric(unlist(dat_outcomes_8to11_f[control_8to11_f[ix_planning_sets_8to11_f],3]))))
plist_8to11_m[[colnum+1]] <- ggplot(df_m, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() + ylim(-2.5,5) + theme_minimal() + theme(legend.position="none") +
  labs(title=paste('Fewer outliers,',outcome_names_8to11[3]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)
p_m <- do.call("grid.arrange", c(plist_8to11_m, nrow=2, top='Ages 8 to 11, Male, Paired Differences'))
p_f <- do.call("grid.arrange", c(plist_8to11_f, nrow=2, top='Ages 8 to 11, Female, Paired Differences'))
ggsave("ages-8to11-male-paireddifs.pdf", p_m, width = 12, height = 9, device = "pdf")
ggsave("ages-8to11-female-paireddifs.pdf", p_f, width = 12, height = 9, device = "pdf")


# Plots of paired differences: Ages 12-17
plist_12to17_m <- plist_12to17_f <- vector(mode="list",length=length(all_outcome_names))
for (colnum in 1:length(all_outcome_names)){
  
  cat(paste('Outcome:',all_outcome_names[colnum]),'\n')
  
  df_m <- melt(data.frame(MatchedDifs=
                            as.numeric(unlist(dat_outcomes_12to17_m[treated_12to17_m[ix_planning_sets_12to17_m],colnum]))-
                            as.numeric(unlist(dat_outcomes_12to17_m[control_12to17_m[ix_planning_sets_12to17_m],colnum]))))
  
  df_f <- melt(data.frame(MatchedDifs=
                            as.numeric(unlist(dat_outcomes_12to17_f[treated_12to17_f[ix_planning_sets_12to17_f],colnum]))-
                            as.numeric(unlist(dat_outcomes_12to17_f[control_12to17_f[ix_planning_sets_12to17_f],colnum]))))
  
  plist_12to17_m[[colnum]] <- ggplot(df_m, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() + theme_minimal() + theme(legend.position="none") +
    labs(title=paste(all_outcome_names[colnum]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)
  
  plist_12to17_f[[colnum]] <- ggplot(df_f, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() + theme_minimal() + theme(legend.position="none") +
    labs(title=paste(all_outcome_names[colnum]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)

  # Two-sided p-values
  cat(paste('Male pvalue:',wilcox.test( df_m$value )$p.value,'\n'))
  cat(paste('Female pvalue:',wilcox.test( df_f$value )$p.value,'\n\n\n'))  
  
}

# Also plot cotinine levels with fewer outliers
df_m <- melt(data.frame(MatchedDifs=
                          as.numeric(unlist(dat_outcomes_12to17_m[treated_12to17_m[ix_planning_sets_12to17_m],5]))-
                          as.numeric(unlist(dat_outcomes_12to17_m[control_12to17_m[ix_planning_sets_12to17_m],5]))))
df_f <- melt(data.frame(MatchedDifs=
                          as.numeric(unlist(dat_outcomes_12to17_f[treated_12to17_f[ix_planning_sets_12to17_f],5]))-
                          as.numeric(unlist(dat_outcomes_12to17_f[control_12to17_f[ix_planning_sets_12to17_f],5]))))
plist_12to17_m[[colnum+1]] <- plist_12to17_m[[6]]
plist_12to17_f[[colnum+1]] <- plist_12to17_f[[6]]

plist_12to17_m[[6]] <- ggplot(df_m, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() + ylim(-10,10) + theme_minimal() + theme(legend.position="none") +
  labs(title=paste('Fewer outliers,',all_outcome_names[5]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)
plist_12to17_f[[6]] <- ggplot(df_f, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() + ylim(-10,10) + theme_minimal() + theme(legend.position="none") +
  labs(title=paste('Fewer outliers,',all_outcome_names[5]),x=NULL,y="")+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)
p_m1 <- do.call("grid.arrange", c(plist_12to17_m[1:6], nrow=2, top='Ages 12 to 17, Male, Paired Differences'))
p_m2 <- do.call("grid.arrange", c(plist_12to17_m[7:11], nrow=2, top='Ages 12 to 17, Male, Paired Differences'))
p_f1 <- do.call("grid.arrange", c(plist_12to17_f[1:6], nrow=2, top='Ages 12 to 17, Female, Paired Differences'))
p_f2 <- do.call("grid.arrange", c(plist_12to17_f[7:11], nrow=2, top='Ages 12 to 17, Female, Paired Differences'))
ggsave("ages-12to17-male-paireddifs1.pdf", p_m1, width = 12, height = 9, device = "pdf")
ggsave("ages-12to17-male-paireddifs2.pdf", p_m2, width = 12, height = 9, device = "pdf")
ggsave("ages-12to17-female-paireddifs1.pdf", p_f1, width = 12, height = 9, device = "pdf")
ggsave("ages-12to17-female-paireddifs2.pdf", p_f2, width = 12, height = 9, device = "pdf")
qpdf::pdf_combine(
  input = c("ages-12to17-male-paireddifs1.pdf","ages-12to17-male-paireddifs2.pdf"),
  output = "ages-12to17-male-paireddifs.pdf"
)
qpdf::pdf_combine(
  input = c("ages-12to17-female-paireddifs1.pdf","ages-12to17-female-paireddifs2.pdf"),
  output = "ages-12to17-female-paireddifs.pdf"
)



################################################################################
#  Estimate optimal contrasts on planning sample
################################################################################

# Find optimal test statistic for each H_k to do testing
mstatlist <- vector('list', length=15)
mstatlist[[1]] <- c(3.5, 0)
mstatlist[[2]] <- c(3.5, 0.5)
mstatlist[[3]] <- c(3.5, 1)
mstatlist[[4]] <- c(3.5, 1.5)
mstatlist[[5]] <- c(3.5, 2)
mstatlist[[6]] <- c(3.5, 2.5)
mstatlist[[7]] <- c(3.5, 3)
mstatlist[[8]] <- c(2.5, 0) # Huber
mstatlist[[9]] <- c(2.5, 0.5) # InnerTrim
mstatlist[[10]] <- c(2.5, 1)
mstatlist[[11]] <- c(2.5, 1.5)
mstatlist[[12]] <- c(2.5, 2)
mstatlist[[13]] <- c(1.5, 0)
mstatlist[[14]] <- c(1.5, 0.5)
mstatlist[[15]] <- c(1.5, 1)

mstat_sensvallist_pilot_8to11_m <- mstat_sensvallist_pilot_8to11_f <- mstat_sensvallist_pilot_8to11 <- 
  mstat_pvallist_pilot_8to11_m <- mstat_pvallist_pilot_8to11_f <- mstat_pvallist_pilot_8to11 <- 
  mstat_sensvallist_pilot_12to17_m <- mstat_sensvallist_pilot_12to17_f <- mstat_sensvallist_pilot_12to17 <- 
  mstat_pvallist_pilot_12to17_m <- mstat_pvallist_pilot_12to17_f <- mstat_pvallist_pilot_12to17 <- 
  vector('list', length=length(mstatlist))
for (mstatix in 1:length(mstatlist)){
  mstat <- mstatlist[[mstatix]]
  
  # Get sensitivity value
  mstat_sensvallist_pilot_8to11_m[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11_m[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11_m[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11_m[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  mstat_sensvallist_pilot_8to11_f[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11_f[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11_f[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11_f[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  mstat_sensvallist_pilot_8to11[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_8to11[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  mstat_sensvallist_pilot_12to17_m[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17_m[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_m[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_m[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_m[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_m[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  mstat_sensvallist_pilot_12to17_f[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17_f[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_f[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_f[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_f[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17_f[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  mstat_sensvallist_pilot_12to17[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval),
    apply(mpdifs_pilot_12to17[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,numGamma=15,alpha=alpha)$sensval))
  
  # Get p-value, used if all Gamma=1
  mstat_pvallist_pilot_8to11_m[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11_m[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11_m[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11_m[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
  
  mstat_pvallist_pilot_8to11_f[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11_f[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11_f[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11_f[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
  
  mstat_pvallist_pilot_8to11[[mstatix]] <- c(
    apply(mpdifs_pilot_8to11[,1:5], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11[,6:14], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_8to11[,15:18], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
  
  mstat_pvallist_pilot_12to17_m[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17_m[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_m[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_m[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_m[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_m[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
  
  mstat_pvallist_pilot_12to17_f[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17_f[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_f[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_f[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_f[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17_f[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
  
  mstat_pvallist_pilot_12to17[[mstatix]] <- c(
    apply(mpdifs_pilot_12to17[,1:2], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17[,3:4], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17[,5:9], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17[,10:18], 2, function(x) main_fixedlam_general(as.matrix(x),-1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval),
    apply(mpdifs_pilot_12to17[,19:22], 2, function(x) main_fixedlam_general(as.matrix(x),1,psi.f=NULL,trim.in.vec=mstat,Gamma=1)$pval))
}


mstat_pilot_8to11_m <- mstat_pilot_8to11_f <- mstat_pilot_8to11 <- vector('list', length=length(mstat_sensvallist_pilot_8to11_m[[1]]))
mstat_pilot_12to17_m <- mstat_pilot_12to17_f <- mstat_pilot_12to17 <- vector('list', length=length(mstat_sensvallist_pilot_12to17_m[[1]]))
for (mstatix in 1:length(mstat_pilot_8to11_m)){
  
  if (max(unlist(lapply(mstat_sensvallist_pilot_8to11_m, function(x) x[mstatix]))) == 1){
    mstat_pilot_8to11_m[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_8to11_m, function(x) x[mstatix]))])
  } else {
    mstat_pilot_8to11_m[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_8to11_m, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_sensvallist_pilot_8to11_f, function(x) x[mstatix]))) == 1){
    mstat_pilot_8to11_f[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_8to11_f, function(x) x[mstatix]))])
  } else {
    mstat_pilot_8to11_f[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_8to11_f, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_sensvallist_pilot_8to11, function(x) x[mstatix]))) == 1){
    mstat_pilot_8to11[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_8to11, function(x) x[mstatix]))])
  } else {
    mstat_pilot_8to11[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_8to11, function(x) x[mstatix]))])
  }

}
for (mstatix in 1:length(mstat_pilot_12to17_m)){

  if (max(unlist(lapply(mstat_sensvallist_pilot_12to17_m, function(x) x[mstatix]))) == 1){
    mstat_pilot_12to17_m[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_12to17_m, function(x) x[mstatix]))])
  } else {
    mstat_pilot_12to17_m[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_12to17_m, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_sensvallist_pilot_12to17_f, function(x) x[mstatix]))) == 1){
    mstat_pilot_12to17_f[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_12to17_f, function(x) x[mstatix]))])
  } else {
    mstat_pilot_12to17_f[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_12to17_f, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_sensvallist_pilot_12to17, function(x) x[mstatix]))) == 1){
    mstat_pilot_12to17[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_pvallist_pilot_12to17, function(x) x[mstatix]))])
  } else {
    mstat_pilot_12to17[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_sensvallist_pilot_12to17, function(x) x[mstatix]))])
  }
  
}

# Note: HEI_total seems less significant them some of its individual components
# Let's make a new, robust dietary composite measure.

# We sometimes have particular outcomes w/ >=50% of paired difs missing; cant get s_k since divides by 0
# We remove those variables k -- if most differences equal zero, there is rarely a difference caused by poverty

# Deal with diet first
# Ages 8-11
mat_pilot_8to11_m_diet <- as.matrix(mpdifs_pilot_8to11_m)[,6:18]
directions_pilot_8to11_m_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_8to11_m_diet <- apply(mat_pilot_8to11_m_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_8to11_m_diet <- mat_pilot_8to11_m_diet[,keep_vars_8to11_m_diet]
directions_pilot_8to11_m_diet <- directions_pilot_8to11_m_diet[keep_vars_8to11_m_diet]
res_pilot_8to11_m_diet <- main_planning_general(mat_pilot_8to11_m_diet,
                                                directions_pilot_8to11_m_diet,
                                                trim.in.list=mstat_pilot_8to11_m[(6:18)[keep_vars_8to11_m_diet]])
mat_pilot_8to11_f_diet <- as.matrix(mpdifs_pilot_8to11_f)[,6:18]
directions_pilot_8to11_f_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_8to11_f_diet <- apply(mat_pilot_8to11_f_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_8to11_f_diet <- mat_pilot_8to11_f_diet[,keep_vars_8to11_f_diet]
directions_pilot_8to11_f_diet <- directions_pilot_8to11_f_diet[keep_vars_8to11_f_diet]
res_pilot_8to11_f_diet <- main_planning_general(mat_pilot_8to11_f_diet,
                                                directions_pilot_8to11_f_diet,
                                                trim.in.list=mstat_pilot_8to11_f[(6:18)[keep_vars_8to11_f_diet]])
mat_pilot_8to11_diet <- as.matrix(mpdifs_pilot_8to11)[,6:18]
directions_pilot_8to11_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_8to11_diet <- apply(mat_pilot_8to11_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_8to11_diet <- mat_pilot_8to11_diet[,keep_vars_8to11_diet]
directions_pilot_8to11_diet <- directions_pilot_8to11_diet[keep_vars_8to11_diet]
res_pilot_8to11_diet <- main_planning_general(mat_pilot_8to11_diet,
                                              directions_pilot_8to11_diet,
                                              trim.in.list=mstat_pilot_8to11[(6:18)[keep_vars_8to11_diet]])
# Ages 12-17
mat_pilot_12to17_m_diet <- as.matrix(mpdifs_pilot_12to17_m)[,10:22]
directions_pilot_12to17_m_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_12to17_m_diet <- apply(mat_pilot_12to17_m_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_12to17_m_diet <- mat_pilot_12to17_m_diet[,keep_vars_12to17_m_diet]
directions_pilot_12to17_m_diet <- directions_pilot_12to17_m_diet[keep_vars_12to17_m_diet]
res_pilot_12to17_m_diet <- main_planning_general(mat_pilot_12to17_m_diet,
                                                 directions_pilot_12to17_m_diet,
                                                 trim.in.list=mstat_pilot_12to17_m[(10:22)[keep_vars_12to17_m_diet]])
mat_pilot_12to17_f_diet <- as.matrix(mpdifs_pilot_12to17_f)[,10:22]
directions_pilot_12to17_f_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_12to17_f_diet <- apply(mat_pilot_12to17_f_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_12to17_f_diet <- mat_pilot_12to17_f_diet[,keep_vars_12to17_f_diet]
directions_pilot_12to17_f_diet <- directions_pilot_12to17_f_diet[keep_vars_12to17_f_diet]
res_pilot_12to17_f_diet <- main_planning_general(mat_pilot_12to17_f_diet,
                                                 directions_pilot_12to17_f_diet,
                                                 trim.in.list=mstat_pilot_12to17_f[(10:22)[keep_vars_12to17_f_diet]])
mat_pilot_12to17_diet <- as.matrix(mpdifs_pilot_12to17)[,10:22]
directions_pilot_12to17_diet <- c(rep('Less',9), rep('Greater',4)) # goes along with listed adequacy and moderation components of HEI
keep_vars_12to17_diet <- apply(mat_pilot_12to17_diet, 2, function(x) mean(x == 0)) < 0.5
mat_pilot_12to17_diet <- mat_pilot_12to17_diet[,keep_vars_12to17_diet]
directions_pilot_12to17_diet <- directions_pilot_12to17_diet[keep_vars_12to17_diet]
res_pilot_12to17_diet <- main_planning_general(mat_pilot_12to17_diet,
                                               directions_pilot_12to17_diet,
                                               trim.in.list=mstat_pilot_12to17[(10:22)[keep_vars_12to17_diet]])


# Combine diet with remaining variables
datalist_pilot_8to11_m <- vector('list',2)
datalist_pilot_8to11_m[[1]]$data <- as.matrix(mpdifs_pilot_8to11_m)[,1:5]
datalist_pilot_8to11_m[[2]]$data <- mat_pilot_8to11_m_diet
datalist_pilot_8to11_m[[2]]$lambda <- res_pilot_8to11_m_diet$lambda
datalist_pilot_8to11_f <- vector('list',2)
datalist_pilot_8to11_f[[1]]$data <- as.matrix(mpdifs_pilot_8to11_f)[,1:5]
datalist_pilot_8to11_f[[2]]$data <- mat_pilot_8to11_f_diet
datalist_pilot_8to11_f[[2]]$lambda <- res_pilot_8to11_f_diet$lambda
datalist_pilot_8to11 <- vector('list',2)
datalist_pilot_8to11[[1]]$data <- as.matrix(mpdifs_pilot_8to11)[,1:5]
datalist_pilot_8to11[[2]]$data <- mat_pilot_8to11_diet
datalist_pilot_8to11[[2]]$lambda <- res_pilot_8to11_diet$lambda

datalist_pilot_12to17_m <- vector('list',2)
datalist_pilot_12to17_m[[1]]$data <- as.matrix(mpdifs_pilot_12to17_m)[,1:9]
datalist_pilot_12to17_m[[2]]$data <- mat_pilot_12to17_m_diet
datalist_pilot_12to17_m[[2]]$lambda <- res_pilot_12to17_m_diet$lambda
datalist_pilot_12to17_f <- vector('list',2)
datalist_pilot_12to17_f[[1]]$data <- as.matrix(mpdifs_pilot_12to17_f)[,1:9]
datalist_pilot_12to17_f[[2]]$data <- mat_pilot_12to17_f_diet
datalist_pilot_12to17_f[[2]]$lambda <- res_pilot_12to17_f_diet$lambda
datalist_pilot_12to17 <- vector('list',2)
datalist_pilot_12to17[[1]]$data <- as.matrix(mpdifs_pilot_12to17)[,1:9]
datalist_pilot_12to17[[2]]$data <- mat_pilot_12to17_diet
datalist_pilot_12to17[[2]]$lambda <- res_pilot_12to17_diet$lambda

# Get optimal statistics to use in combination
mstat_diet_sensvallist_pilot_8to11_m <- mstat_diet_sensvallist_pilot_8to11_f <- mstat_diet_sensvallist_pilot_8to11 <- 
  mstat_diet_pvallist_pilot_8to11_m <- mstat_diet_pvallist_pilot_8to11_f <- mstat_diet_pvallist_pilot_8to11 <- 
  mstat_diet_sensvallist_pilot_12to17_m <- mstat_diet_sensvallist_pilot_12to17_f <- mstat_diet_sensvallist_pilot_12to17 <- 
  mstat_diet_pvallist_pilot_12to17_m <- mstat_diet_pvallist_pilot_12to17_f <- mstat_diet_pvallist_pilot_12to17 <- 
  vector('list', length=length(mstatlist))
for (mstatix in 1:length(mstatlist)){
  mstat <- mstatlist[[mstatix]]
  
  # Get sensitivity value
  mstat_diet_sensvallist_pilot_8to11_m[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11_m[2],
                                                                           1,
                                                                           psi.f=NULL,
                                                                           trim.in.vec=mstat,
                                                                           numGamma=15,
                                                                           alpha=alpha)$sensval
  mstat_diet_sensvallist_pilot_8to11_f[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11_f[2],
                                                                           1,
                                                                           psi.f=NULL,
                                                                           trim.in.vec=mstat,
                                                                           numGamma=15,
                                                                           alpha=alpha)$sensval  
  mstat_diet_sensvallist_pilot_8to11[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11[2],
                                                                         1,
                                                                         psi.f=NULL,
                                                                         trim.in.vec=mstat,
                                                                         numGamma=15,
                                                                         alpha=alpha)$sensval  
  mstat_diet_sensvallist_pilot_12to17_m[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17_m[2],
                                                                            1,
                                                                            psi.f=NULL,
                                                                            trim.in.vec=mstat,
                                                                            numGamma=15,
                                                                            alpha=alpha)$sensval  
  mstat_diet_sensvallist_pilot_12to17_f[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17_f[2],
                                                                            1,
                                                                            psi.f=NULL,
                                                                            trim.in.vec=mstat,
                                                                            numGamma=15,
                                                                            alpha=alpha)$sensval  
  mstat_diet_sensvallist_pilot_12to17[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17[2],
                                                                          1,
                                                                          psi.f=NULL,
                                                                          trim.in.vec=mstat,
                                                                          numGamma=15,
                                                                          alpha=alpha)$sensval
  # Get p-value, used if all Gamma=1
  mstat_diet_pvallist_pilot_8to11_m[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11_m[2],
                                                                           1,
                                                                           psi.f=NULL,
                                                                           trim.in.vec=mstat,
                                                                           Gamma=1)$pval
  mstat_diet_pvallist_pilot_8to11_f[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11_f[2],
                                                                           1,
                                                                           psi.f=NULL,
                                                                           trim.in.vec=mstat,
                                                                           Gamma=1)$pval  
  mstat_diet_pvallist_pilot_8to11[[mstatix]] <- main_fixedlam_general(datalist_pilot_8to11[2],
                                                                         1,
                                                                         psi.f=NULL,
                                                                         trim.in.vec=mstat,
                                                                         Gamma=1)$pval  
  mstat_diet_pvallist_pilot_12to17_m[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17_m[2],
                                                                            1,
                                                                            psi.f=NULL,
                                                                            trim.in.vec=mstat,
                                                                            Gamma=1)$pval  
  mstat_diet_pvallist_pilot_12to17_f[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17_f[2],
                                                                            1,
                                                                            psi.f=NULL,
                                                                            trim.in.vec=mstat,
                                                                            Gamma=1)$pval  
  mstat_diet_pvallist_pilot_12to17[[mstatix]] <- main_fixedlam_general(datalist_pilot_12to17[2],
                                                                          1,
                                                                          psi.f=NULL,
                                                                          trim.in.vec=mstat,
                                                                          Gamma=1)$pval
}
mstat_diet_pilot_8to11_m <- mstat_diet_pilot_8to11_f <- mstat_diet_pilot_8to11 <- vector('list', length=length(mstat_diet_sensvallist_pilot_8to11_m[[1]]))
mstat_diet_pilot_12to17_m <- mstat_diet_pilot_12to17_f <- mstat_diet_pilot_12to17 <- vector('list', length=length(mstat_diet_sensvallist_pilot_12to17_m[[1]]))
for (mstatix in 1:length(mstat_diet_pilot_8to11_m)){
  
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_8to11_m, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_8to11_m[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_8to11_m, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_8to11_m[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_8to11_m, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_8to11_f, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_8to11_f[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_8to11_f, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_8to11_f[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_8to11_f, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_8to11, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_8to11[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_8to11, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_8to11[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_8to11, function(x) x[mstatix]))])
  }
  
}
for (mstatix in 1:length(mstat_diet_pilot_12to17_m)){
  
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_12to17_m, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_12to17_m[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_12to17_m, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_12to17_m[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_12to17_m, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_12to17_f, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_12to17_f[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_12to17_f, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_12to17_f[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_12to17_f, function(x) x[mstatix]))])
  }
  if (max(unlist(lapply(mstat_diet_sensvallist_pilot_12to17, function(x) x[mstatix]))) == 1){
    mstat_diet_pilot_12to17[[mstatix]] <- 
      unlist(mstatlist[which.min(lapply(mstat_diet_pvallist_pilot_12to17, function(x) x[mstatix]))])
  } else {
    mstat_diet_pilot_12to17[[mstatix]] <- 
      unlist(mstatlist[which.max(lapply(mstat_diet_sensvallist_pilot_12to17, function(x) x[mstatix]))])
  }
  
}

# Do estimation
res_pilot_8to11_m <- main_planning_general(datalist_pilot_8to11_m,
                                           rep('Greater',6),
                                           trim.in.list=list(mstat_pilot_8to11_m[(1:5)], mstat_diet_pilot_8to11_m))
res_pilot_8to11_f <- main_planning_general(datalist_pilot_8to11_f,
                                           rep('Greater',6),
                                           trim.in.list=list(mstat_pilot_8to11_f[(1:5)], mstat_diet_pilot_8to11_f))
res_pilot_8to11 <- main_planning_general(datalist_pilot_8to11,
                                         rep('Greater',6),
                                         trim.in.list=list(mstat_pilot_8to11[(1:5)], mstat_diet_pilot_8to11))
res_pilot_12to17_m <- main_planning_general(datalist_pilot_12to17_m,
                                            c(rep('Greater',2), rep('Less',2), rep('Greater',6)), # worry about lower phys activity
                                            trim.in.list=list(mstat_pilot_12to17_m[(1:9)], mstat_diet_pilot_12to17_m))
res_pilot_12to17_f <- main_planning_general(datalist_pilot_12to17_f,
                                            c(rep('Greater',2), rep('Less',2), rep('Greater',6)), # worry about lower phys activity
                                            trim.in.list=list(mstat_pilot_12to17_f[(1:9)], mstat_diet_pilot_12to17_f))
res_pilot_12to17 <- main_planning_general(datalist_pilot_12to17,
                                          c(rep('Greater',2), rep('Less',2), rep('Greater',6)), # worry about lower phys activity
                                          trim.in.list=list(mstat_pilot_12to17[(1:9)], mstat_diet_pilot_12to17))



################################################################################
#  Do testing on analysis sample
################################################################################
# Get paired differences on analysis sample
mpdifs_analysis_8to11_m <- dat_outcomes_8to11_m[treated_8to11_m[ix_analysis_sets_8to11_m],]-dat_outcomes_8to11_m[control_8to11_m[ix_analysis_sets_8to11_m],]
mpdifs_analysis_8to11_f <- dat_outcomes_8to11_f[treated_8to11_f[ix_analysis_sets_8to11_f],]-dat_outcomes_8to11_f[control_8to11_f[ix_analysis_sets_8to11_f],]
mpdifs_analysis_12to17_m <- dat_outcomes_12to17_m[treated_12to17_m[ix_analysis_sets_12to17_m],]-dat_outcomes_12to17_m[control_12to17_m[ix_analysis_sets_12to17_m],]
mpdifs_analysis_12to17_f <- dat_outcomes_12to17_f[treated_12to17_f[ix_analysis_sets_12to17_f],]-dat_outcomes_12to17_f[control_12to17_f[ix_analysis_sets_12to17_f],]
mpdifs_analysis_8to11 <- rbind(mpdifs_analysis_8to11_m, mpdifs_analysis_8to11_f)
mpdifs_analysis_12to17 <- rbind(mpdifs_analysis_12to17_m, mpdifs_analysis_12to17_f)




gamma_vector <- c(1, 1.05, 1.15, 1.25, 1.5, 1.75, 2, 2.25, 2.5)

for (gamma in gamma_vector){
  
  # Get p-values for each H_k (i.e., leaf node p-values) using one-sided tests
  ind_pvals_list <- vector("list", length=4)
  names(ind_pvals_list) <- c('Age 8-11, Male',
                             'Age 8-11, Female',
                             'Age 12-17, Male',
                             'Age 12-17, Female')
  # Ages 8-11, Male
  # Calculate diet p-value, include remaining variables
  mat_analysis_8to11_m_diet <- as.matrix(mpdifs_analysis_8to11_m[,6:18])
  mat_analysis_8to11_m_diet <- mat_analysis_8to11_m_diet[,keep_vars_8to11_m_diet]
  indp_analysis_8to11_m_diet <- main_fixedlam_general(mat_analysis_8to11_m_diet,
                                                      res_pilot_8to11_m_diet$lambda,
                                                      trim.in.vec=unlist(mstat_diet_pilot_8to11_m),
                                                      Gamma=gamma)$pval
  ind_pvals_list[[1]] <- c(mapply(
    FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
    col = as.data.frame(mpdifs_analysis_8to11_m[,1:5]),
    trim_param = mstat_pilot_8to11_m[1:5],
    SIMPLIFY = F
    ), indp_analysis_8to11_m_diet)
  names(ind_pvals_list[[1]])[6] <- 'Diet'
  # Ages 8-11, Female
  # Calculate diet p-value, include remaining variables
  mat_analysis_8to11_f_diet <- as.matrix(mpdifs_analysis_8to11_f[,6:18])
  mat_analysis_8to11_f_diet <- mat_analysis_8to11_f_diet[,keep_vars_8to11_f_diet]
  indp_analysis_8to11_f_diet <- main_fixedlam_general(mat_analysis_8to11_f_diet,
                                                      res_pilot_8to11_f_diet$lambda,
                                                      trim.in.vec=unlist(mstat_diet_pilot_8to11_f),
                                                      Gamma=gamma)$pval
  ind_pvals_list[[2]] <- c(mapply(
    FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
    col = as.data.frame(mpdifs_analysis_8to11_f[,1:5]),
    trim_param = mstat_pilot_8to11_f[1:5],
    SIMPLIFY = F
  ), indp_analysis_8to11_f_diet)
  names(ind_pvals_list[[2]])[6] <- 'Diet'
  # Ages 12-17, Male
  # Calculate diet p-value, include remaining variables
  mat_analysis_12to17_m_diet <- as.matrix(mpdifs_analysis_12to17_m[,10:22])
  mat_analysis_12to17_m_diet <- mat_analysis_12to17_m_diet[,keep_vars_12to17_m_diet]
  indp_analysis_12to17_m_diet <- main_fixedlam_general(mat_analysis_12to17_m_diet,
                                                       res_pilot_12to17_m_diet$lambda,
                                                       trim.in.vec=unlist(mstat_diet_pilot_12to17_m),
                                                       Gamma=gamma)$pval
  ind_pvals_list[[3]] <- 
    c(mapply(
    FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
    col = as.data.frame(mpdifs_analysis_12to17_m[,1:2]),
    trim_param = mstat_pilot_12to17_m[1:2],
    SIMPLIFY = F),
    mapply(
      FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),-1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
      col = as.data.frame(mpdifs_analysis_12to17_m[,3:4]),
      trim_param = mstat_pilot_12to17_m[3:4],
      SIMPLIFY = F),
    mapply(
      FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
      col = as.data.frame(mpdifs_analysis_12to17_m[,5:9]),
      trim_param = mstat_pilot_12to17_m[5:9],
      SIMPLIFY = F),
    indp_analysis_12to17_m_diet)
  names(ind_pvals_list[[3]])[10] <- 'Diet'
  # Ages 12-17, Female
  # Calculate diet p-value, include remaining variables
  mat_analysis_12to17_f_diet <- as.matrix(mpdifs_analysis_12to17_f[,10:22])
  mat_analysis_12to17_f_diet <- mat_analysis_12to17_f_diet[,keep_vars_12to17_f_diet]
  indp_analysis_12to17_f_diet <- main_fixedlam_general(mat_analysis_12to17_f_diet,
                                                       res_pilot_12to17_f_diet$lambda,
                                                       trim.in.vec=unlist(mstat_diet_pilot_12to17_f),
                                                       Gamma=gamma)$pval
  ind_pvals_list[[4]] <- 
    c(mapply(
      FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
      col = as.data.frame(mpdifs_analysis_12to17_f[,1:2]),
      trim_param = mstat_pilot_12to17_f[1:2],
      SIMPLIFY = F),
      mapply(
        FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),-1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
        col = as.data.frame(mpdifs_analysis_12to17_f[,3:4]),
        trim_param = mstat_pilot_12to17_f[3:4],
        SIMPLIFY = F),
      mapply(
        FUN = function(col, trim_param) {main_fixedlam_general(as.matrix(col),1,trim.in.list=list(trim_param),Gamma=gamma)$pval},
        col = as.data.frame(mpdifs_analysis_12to17_f[,5:9]),
        trim_param = mstat_pilot_12to17_f[5:9],
        SIMPLIFY = F),
      indp_analysis_12to17_f_diet)
  names(ind_pvals_list[[4]])[10] <- 'Diet'
  
  # Get global null p-values using directions we estimated on pilot data
  # Ages 8-11
  datalist_analysis_8to11_m <- vector('list',2)
  datalist_analysis_8to11_m[[1]]$data <- as.matrix(mpdifs_analysis_8to11_m)[,1:5]
  datalist_analysis_8to11_m[[2]]$data <- mat_analysis_8to11_m_diet
  datalist_analysis_8to11_m[[2]]$lambda <- res_pilot_8to11_m_diet$lambda
  res_analysis_8to11_m <- main_fixedlam_general(datalist_analysis_8to11_m,
                                                res_pilot_8to11_m$lambda,
                                                trim.in.list=list(mstat_pilot_8to11_m[(1:5)], mstat_diet_pilot_8to11_m),
                                                Gamma=gamma)$pval
  datalist_analysis_8to11_f <- vector('list',2)
  datalist_analysis_8to11_f[[1]]$data <- as.matrix(mpdifs_analysis_8to11_f)[,1:5]
  datalist_analysis_8to11_f[[2]]$data <- mat_analysis_8to11_f_diet
  datalist_analysis_8to11_f[[2]]$lambda <- res_pilot_8to11_f_diet$lambda
  res_analysis_8to11_f <- main_fixedlam_general(datalist_analysis_8to11_f,
                                                res_pilot_8to11_f$lambda,
                                                trim.in.list=list(mstat_pilot_8to11_f[(1:5)], mstat_diet_pilot_8to11_f),
                                                Gamma=gamma)$pval
  mat_analysis_8to11_diet <- as.matrix(mpdifs_analysis_8to11[,6:18])
  mat_analysis_8to11_diet <- mat_analysis_8to11_diet[,keep_vars_8to11_diet]
  datalist_analysis_8to11 <- vector('list',2)
  datalist_analysis_8to11[[1]]$data <- as.matrix(mpdifs_analysis_8to11)[,1:5]
  datalist_analysis_8to11[[2]]$data <- mat_analysis_8to11_diet
  datalist_analysis_8to11[[2]]$lambda <- res_pilot_8to11_diet$lambda
  res_analysis_8to11 <- main_fixedlam_general(datalist_analysis_8to11,
                                              res_pilot_8to11$lambda,
                                              trim.in.list=list(mstat_pilot_8to11[(1:5)], mstat_diet_pilot_8to11),
                                              Gamma=gamma)$pval
  # Ages 12-17
  datalist_analysis_12to17_m <- vector('list',2)
  datalist_analysis_12to17_m[[1]]$data <- as.matrix(mpdifs_analysis_12to17_m)[,1:9]
  datalist_analysis_12to17_m[[2]]$data <- mat_analysis_12to17_m_diet
  datalist_analysis_12to17_m[[2]]$lambda <- res_pilot_12to17_m_diet$lambda
  res_analysis_12to17_m <- main_fixedlam_general(datalist_analysis_12to17_m,
                                                 res_pilot_12to17_m$lambda,
                                                 trim.in.list=list(mstat_pilot_12to17_m[(1:9)], mstat_diet_pilot_12to17_m),
                                                 Gamma=gamma)$pval
  datalist_analysis_12to17_f <- vector('list',2)
  datalist_analysis_12to17_f[[1]]$data <- as.matrix(mpdifs_analysis_12to17_f)[,1:9]
  datalist_analysis_12to17_f[[2]]$data <- mat_analysis_12to17_f_diet
  datalist_analysis_12to17_f[[2]]$lambda <- res_pilot_12to17_f_diet$lambda
  res_analysis_12to17_f <- main_fixedlam_general(datalist_analysis_12to17_f,
                                                 res_pilot_12to17_f$lambda,
                                                 trim.in.list=list(mstat_pilot_12to17_f[(1:9)], mstat_diet_pilot_12to17_f),
                                                 Gamma=gamma)$pval
  mat_analysis_12to17_diet <- as.matrix(mpdifs_analysis_12to17[,10:22])
  mat_analysis_12to17_diet <- mat_analysis_12to17_diet[,keep_vars_12to17_diet]
  datalist_analysis_12to17 <- vector('list',2)
  datalist_analysis_12to17[[1]]$data <- as.matrix(mpdifs_analysis_12to17)[,1:9]
  datalist_analysis_12to17[[2]]$data <- mat_analysis_12to17_diet
  datalist_analysis_12to17[[2]]$lambda <- res_pilot_12to17_diet$lambda
  res_analysis_12to17 <- main_fixedlam_general(datalist_analysis_12to17,
                                               res_pilot_12to17$lambda,
                                               trim.in.list=list(mstat_pilot_12to17[(1:9)], mstat_diet_pilot_12to17),
                                               Gamma=gamma)$pval
  
  # Get root p-value
  if (any(c(res_analysis_8to11,
            res_analysis_12to17)==0)
  ){
    root_pval=0
  } else { # run Fisher combination on p-values to get root value
    root_pval <- truncatedP(c(res_analysis_8to11,
                              res_analysis_12to17),
                            trunc=1)
  }
  
  # Aggregate all p-values and run the inheritance procedure (Goeman & Finos, 2012)
  nodes <- list(
    #### ROOT NODE
    Root = list(parents = character(0),
                children = c("Ages8to11","Ages12to17"),
                p = root_pval),
    
    #### AGEGRP NODES
    Ages8to11 = list(parents = "Root",
                     children = c("Ages8to11_M","Ages8to11_F"),
                     p = res_analysis_8to11),
    Ages12to17 = list(parents = "Root",
                      children = c("Ages12to17_M","Ages12to17_F"),
                      p = res_analysis_12to17),
    
    #### AGEGRP/SEX NODES
    Ages8to11_M = list(parents = "Ages8to11",
                       children = paste0("Ages8to11_M","_",c(outcome_names_8to11[1:5],'Diet')),
                       p = res_analysis_8to11_m),
    Ages8to11_F = list(parents = "Ages8to11",
                       children = paste0("Ages8to11_F","_",c(outcome_names_8to11[1:5],'Diet')),
                       p = res_analysis_8to11_f),
    Ages12to17_M = list(parents = "Ages12to17",
                        children = paste0("Ages12to17_M","_",c(all_outcome_names[1:9],'Diet')),
                        p = res_analysis_12to17_m),
    Ages12to17_F = list(parents = "Ages12to17",
                        children = paste0("Ages12to17_F","_",c(all_outcome_names[1:9],'Diet')),
                        p = res_analysis_12to17_f),
    
    #### OUTCOME NODES (LEAVES)                  
    Ages8to11_M_BMXBMI = list(parents = "Ages8to11_M",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['BMXBMI'] )),
    Ages8to11_M_waist_height = list(parents = "Ages8to11_M",
                                    children = character(0),
                                    p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['waist_height'] )),
    Ages8to11_M_LBXCOT = list(parents = "Ages8to11_M",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['LBXCOT'] )),
    Ages8to11_M_sys_bp = list(parents = "Ages8to11_M",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['sys_bp'] )),
    Ages8to11_M_non_HDL_chol = list(parents = "Ages8to11_M",
                                    children = character(0),
                                    p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['non_HDL_chol'] )),
    Ages8to11_M_Diet = list(parents = "Ages8to11_M",
                            children = character(0),
                            p = as.numeric( ind_pvals_list[['Age 8-11, Male']]['Diet'] )),
    
    Ages8to11_F_BMXBMI = list(parents = "Ages8to11_F",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['BMXBMI'] )),
    Ages8to11_F_waist_height = list(parents = "Ages8to11_F",
                                    children = character(0),
                                    p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['waist_height'] )),
    Ages8to11_F_LBXCOT = list(parents = "Ages8to11_F",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['LBXCOT'] )),
    Ages8to11_F_sys_bp = list(parents = "Ages8to11_F",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['sys_bp'] )),
    Ages8to11_F_non_HDL_chol = list(parents = "Ages8to11_F",
                                    children = character(0),
                                    p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['non_HDL_chol'] )),
    Ages8to11_F_Diet = list(parents = "Ages8to11_F",
                            children = character(0),
                            p = as.numeric( ind_pvals_list[['Age 8-11, Female']]['Diet'] )),
    
    Ages12to17_M_BMXBMI = list(parents = "Ages12to17_M",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['BMXBMI'] )),
    Ages12to17_M_waist_height = list(parents = "Ages12to17_M",
                                     children = character(0),
                                     p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['waist_height'] )),
    Ages12to17_M_moderate_leisure = list(parents = "Ages12to17_M",
                                         children = character(0),
                                         p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['moderate_leisure'] )),
    Ages12to17_M_vigorous_leisure = list(parents = "Ages12to17_M",
                                         children = character(0),
                                         p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['vigorous_leisure'] )),
    Ages12to17_M_LBXCOT = list(parents = "Ages12to17_M",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['LBXCOT'] )),
    Ages12to17_M_sys_bp = list(parents = "Ages12to17_M",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['sys_bp'] )),
    Ages12to17_M_LBXGH = list(parents = "Ages12to17_M",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['LBXGH'] )),
    Ages12to17_M_non_HDL_chol = list(parents = "Ages12to17_M",
                                     children = character(0),
                                     p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['non_HDL_chol'] )),
    Ages12to17_M_eGFR = list(parents = "Ages12to17_M",
                             children = character(0),
                             p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['eGFR'] )),
    Ages12to17_M_Diet = list(parents = "Ages12to17_M",
                             children = character(0),
                             p = as.numeric( ind_pvals_list[['Age 12-17, Male']]['Diet'] )),
    
    Ages12to17_F_BMXBMI = list(parents = "Ages12to17_F",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['BMXBMI'] )),
    Ages12to17_F_waist_height = list(parents = "Ages12to17_F",
                                     children = character(0),
                                     p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['waist_height'] )),
    Ages12to17_F_moderate_leisure = list(parents = "Ages12to17_F",
                                         children = character(0),
                                         p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['moderate_leisure'] )),
    Ages12to17_F_vigorous_leisure = list(parents = "Ages12to17_F",
                                         children = character(0),
                                         p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['vigorous_leisure'] )),
    Ages12to17_F_LBXCOT = list(parents = "Ages12to17_F",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['LBXCOT'] )),
    Ages12to17_F_sys_bp = list(parents = "Ages12to17_F",
                               children = character(0),
                               p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['sys_bp'] )),
    Ages12to17_F_LBXGH = list(parents = "Ages12to17_F",
                              children = character(0),
                              p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['LBXGH'] )),
    Ages12to17_F_non_HDL_chol = list(parents = "Ages12to17_F",
                                     children = character(0),
                                     p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['non_HDL_chol'] )),
    Ages12to17_F_eGFR = list(parents = "Ages12to17_F",
                             children = character(0),
                             p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['eGFR'] )),
    Ages12to17_F_Diet = list(parents = "Ages12to17_F",
                             children = character(0),
                             p = as.numeric( ind_pvals_list[['Age 12-17, Female']]['Diet'] ))
  )
  # Do not test for nodes k st lambda_k = 0 on pilot sample
  prune_sex_node <- function(nodes, sex_node, lambda_vec) {
    children <- nodes[[sex_node]]$children
    if (length(children) != length(lambda_vec)) {
      stop(paste0("Lambda length mismatch in sex node ", sex_node,
                  ": length(children) = ", length(children),
                  ", length(lambda) = ", length(lambda_vec)))
    }
    # which indices to drop
    drop_idx <- which(round(lambda_vec, 16) == 0)
    if (length(drop_idx) == 0) return(nodes)
    # actual node names to remove
    drop_nodes <- children[drop_idx]
    # remove nodes from the list
    nodes <- nodes[ setdiff(names(nodes), drop_nodes) ]
    # update children list of this sex node
    nodes[[sex_node]]$children <- setdiff(children, drop_nodes)
    return(nodes)
  }
  # Prune Ages 8–11
  nodes <- prune_sex_node(
    nodes,
    sex_node = "Ages8to11_M",
    lambda_vec = res_pilot_8to11_m$lambda
  )
  nodes <- prune_sex_node(
    nodes,
    sex_node = "Ages8to11_F",
    lambda_vec = res_pilot_8to11_f$lambda
  )
  # Prune Ages 12–17
  nodes <- prune_sex_node(
    nodes,
    sex_node = "Ages12to17_M",
    lambda_vec = res_pilot_12to17_m$lambda
  )
  nodes <- prune_sex_node(
    nodes,
    sex_node = "Ages12to17_F",
    lambda_vec = res_pilot_12to17_f$lambda
  )
  
  nodes <- nodes[as.logical(unlist(lapply(nodes,function(x)!is.null(x))))]
  
  # Un-adjusted p-values
  write_json(
    list(list(
      gamma = gamma,
      pvalues = as.list(sort(unlist(lapply(nodes, function(x) x$p))))
    )),
    paste0("r_output/unadj_pvalues_gamma",sprintf('%1.2f',gamma),".json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  
  # Get adjusted p-values
  adjp <- run_inheritance(nodes, Shaffer = TRUE, homogeneous = FALSE)
  write_json(
    list(list(
      gamma = gamma,
      pvalues = as.list(sort(adjp))
    )),
    paste0("r_output/pvalues_gamma",sprintf('%1.2f',gamma),".json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  sort(round(adjp,3))


}






qpdf::pdf_combine(
  input = paste0('outcomes_tree_gamma',sprintf('%1.2f',gamma_vector),'.pdf'),
  output = "outcomes_tree_gamma_combined.pdf"
)




# # Examine uncorrected pvals without correction for multiple testing
# lapply(ind_pvals_list, function(x) round(x,3))
# res_analysis_8to11_m
# res_analysis_8to11_f
# res_analysis_8to11
# res_analysis_12to17_m
# res_analysis_12to17_f
# res_analysis_12to17
# root_pval
# 
# 
# .05 / (length(adjp) - 7) # for regular Bonf. correction without inheritance
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### full sample differences
# mpdifs_full_8to11_m <- dat_outcomes_8to11_m[treated_8to11_m,]-dat_outcomes_8to11_m[control_8to11_m,]
# mpdifs_full_8to11_f <- dat_outcomes_8to11_f[treated_8to11_f,]-dat_outcomes_8to11_f[control_8to11_f,]
# mpdifs_full_12to17_m <- dat_outcomes_12to17_m[treated_12to17_m,]-dat_outcomes_12to17_m[control_12to17_m,]
# mpdifs_full_12to17_f <- dat_outcomes_12to17_f[treated_12to17_f,]-dat_outcomes_12to17_f[control_12to17_f,]
# mpdifs_full_8to11 <- rbind(mpdifs_full_8to11_m, mpdifs_full_8to11_f)
# mpdifs_full_12to17 <- rbind(mpdifs_full_12to17_m, mpdifs_full_12to17_f)
# 
# 
# 
# 
# 
# 
