p.rows = c(2,5,25,100,150)
normvec = chisqvec = chibarsqvec = corr.chibarsqvec = numeric()
alpha = 0.05

for (p in p.rows){
  normvec <- c(normvec, qnorm(1-alpha))
  chisqvec <- c(chisqvec, sqrt(qchisq(1-alpha, p)))
  chibarsqvec <- c(chibarsqvec, sqrt(qchibarsq(1-alpha, diag(p)))) # assumes corr=0
}

normvec
chisqvec
chibarsqvec

### use fast-weights.R for quantiles of chibarsq dist. w/ nonzero correlation
