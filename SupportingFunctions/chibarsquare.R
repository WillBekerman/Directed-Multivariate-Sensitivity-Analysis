################################################################################
#       Cumulative Distribution Function for Chi-Bar-Squared distribution
################################################################################

pchibarsq=function(q, V, lower.tail=TRUE, log.p=FALSE)
{
  n<-nrow(V)
  
  wts<-wchibarsq(V) 
  
  ans<-pchisq(q, 0, lower.tail=FALSE)*wts[1L]+
    pchisq(q, n, lower.tail=FALSE)*wts[n+1L]
  for(i in seq_len(n-1)){
    ans<-ans+pchisq(q, i, lower.tail=FALSE)*wts[i+1L]
  }
  ans[q<=0]<-1
  
  ans <- if(isTRUE(lower.tail)) 1-ans else ans
  if(isTRUE(log.p)) log(ans) else ans
}

################################################################################
#       Weights Function for Chi-Bar-Squared distribution
################################################################################

wchibarsq=function(V) ## weights
{## Ref: Kudo A. 1963. Biometrika, 50, pg 403  (page 414)
  stopifnot(is.matrix(V) && diff(dim(V))==0L && nrow(V)>0L)
  n<-nrow(V)
  P<-function(idx){## V, n, i
    pmvnorm(rep(0,n-i), rep(Inf, n-i), sigma=solve(V[-idx,-idx,drop=FALSE]))*
      pmvnorm(rep(0, i), rep(Inf, i), 
              sigma=V[idx,idx,drop=FALSE]-V[idx, -idx, drop=FALSE]%*%solve(V[-idx,-idx,drop=FALSE], V[-idx, idx, drop=FALSE])
      )
  }
  ans<-numeric(n+1L)
  ans[1]<-pmvnorm(rep(0,n),rep(Inf,n),sigma=solve(V))[[1]]
  ans[n+1L]<-pmvnorm(rep(0,n), rep(Inf,n), sigma=V)[[1]]
  for(i in seq(1L, n-1L, by=1L)) ans[i+1]=sum(combn(n, i, P))
  ans
}

################################################################################
#       Quantile Function for Chi-Bar-Squared distribution
################################################################################
qchibarsq =function(p, V)
{
  K = ncol(V)
  pchiroot = function(q, V){pchibarsq(q, V) - p}
  uniroot(pchiroot, c(0, qchisq(p, K)), V=V)$root
  
}

