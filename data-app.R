setwd("C:/Users/bekerman/Downloads/relevant-xpt-files")

library(tidyverse)
library(MASS)
library(haven)
library(caret)
library(optmatch)

#' Runs propensity score matching on data.
#'
#' @param data_matching .
#' @param data_outcomes .
#' @param treated Vector of treated indices.
#' @param calipersd Propensity score caliper.
#' @return .
run_PSM <- function(data_matching,data_outcomes,treated,calipersd){
  ## Function for computing rank based Mahalanobis distance.
  smahal=function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-as.matrix(X[z==0,])
    Xt<-as.matrix(X[z==1,])
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }
  ## Function for adding a propensity score caliper to a distance matrix dmat
  addcaliper=function(dmat,z,logitp,calipersd,penalty=1000){
    # Pooled within group standard devation
    sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
    adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
    adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
    dmat=dmat+adif*penalty
    dmat
  }
  datatemp=as.data.frame(cbind(data_matching,treated))
  names(datatemp)[ncol(datatemp)] <- 'treated'
  ## Propensity score model
  propscore.model=glm(treated ~ ., data = datatemp, family = "binomial")
  # pscoredat=datatemp
  # #split training (80%) and testing (20%)
  # parts = createDataPartition(pscoredat$treated,p=.80,list=F)
  # train = pscoredat[parts,]; test = pscoredat[-parts,]
  # #define predictor and response variables in training set
  # train_x = data.matrix(subset(train,select=-c(treated)));train_y = train$treated
  # #define predictor and response variables in testing set
  # test_x = data.matrix(subset(test,select=-c(treated))); test_y = test$treated
  # #define final training and testing sets
  # xgb_train = xgb.DMatrix(data=train_x,label=train_y); xgb_test = xgb.DMatrix(data=test_x,label=test_y)
  # #define watchlist
  # watchlist = list(train=xgb_train, test=xgb_test)
  # #fit XGBoost model and display training and testing data at each round
  # model = xgb.train(data = xgb_train, max.depth = 3, watchlist=watchlist, objective='binary:logistic', nrounds = 10, verbose = 0)
  # #define final model
  # propscore.model = xgboost(data = xgb_train, max.depth = 3, objective='binary:logistic', nrounds = which.min(model$evaluation_log$test_logloss), verbose = 0)
  dmy=dummyVars('treated~.',data=datatemp)
  Xmat=data.frame(predict(dmy,newdata=datatemp))
  # Xmat=subset(datatemp,select=-c(treated))
  Xmatmahal=Xmat
  datatemp$logit.ps=as.numeric(predict(propscore.model,newx=as.matrix(Xmat)))
  # datatemp$logit.ps=as.numeric(predict(propscore.model,newdata=xgb.DMatrix(sapply(Xmat,as.numeric))))
  ## Use Hansen (2009)’s rule for removing subjects who lack overlap
  logit.propscore=datatemp$logit.ps
  pooled.sd.logit.propscore=sqrt(var(logit.propscore[datatemp$treated==1])/2+var
                                 (logit.propscore[datatemp$treated==0])/2)
  min.treated.logit.propscore=min(logit.propscore[datatemp$treated==1])
  max.control.logit.propscore=max(logit.propscore[datatemp$treated==0])
  ## How many treated and control subjects lack overlap by Hansen's criterion
  no.treated.lack.overlap=0#sum(logit.propscore[datatemp$treated==1]>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))
  no.control.lack.overlap=0#sum(logit.propscore[datatemp$treated==0]<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore))
  ## If there are subjects who lack overlap, remove them from the datatemp dataset
  datatemp.original=datatemp
  datatemp.full=datatemp
  Xmat.original=Xmat
  Xmat.full=Xmat
  which.remove=NA
  # if(no.treated.lack.overlap+no.control.lack.overlap>0){
  #   which.remove=which((logit.propscore>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))|(logit.propscore<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore)))
  #   datatemp=datatemp[-which.remove,]
  #   datatemp.full=rbind(datatemp,datatemp.original[which.remove,])
  #   Xmat=Xmat[-which.remove,]
  #   Xmat.full=rbind(Xmat,Xmat.original[which.remove,])
  #   Xmatmahal=Xmatmahal[-which.remove,]
  #   logit.propscore=c(logit.propscore[-which.remove],logit.propscore[which.remove])
  # }
  ## For the purposes of balance checking later, in datatemp.full, append
  ## the removed rows of datatemp to the end of datatemp
  rownames(datatemp)=seq(1,nrow(datatemp),1)
  Xmatmahal$logit.ps=datatemp$logit.ps
  ## Rank based Mahalanobis distance with caliper
  distmat=smahal(datatemp$treated,as.matrix(Xmatmahal$logit.ps))
  distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=calipersd)
  rownames(distmat)=rownames(datatemp)[datatemp$treated==1]
  colnames(distmat)=rownames(datatemp)[datatemp$treated==0]
  ## Match nocontrols.per.match to each treated unit
  nocontrols.per.match=1
  matchvec=pairmatch(distmat,controls=nocontrols.per.match,data=datatemp)
  datatemp$matchvec=matchvec
  matchedset.index=substr(matchvec,start=3,stop=10)
  matchedset.index.numeric=as.numeric(matchedset.index)
  ## Have matchedset.index.numeric.full append 0 to matchedset.index.numeric for
  ## the removed subjects
  if(no.control.lack.overlap+no.treated.lack.overlap==0){
    matchedset.index.numeric.full=matchedset.index.numeric
  }
  if(no.control.lack.overlap+no.treated.lack.overlap>0){
    matchedset.index.numeric.full=c(matchedset.index.numeric,rep(0,no.control.lack.overlap+no.treated.lack.overlap))
  }
  ## Create a matrix saying which control units each treated unit is matched to
  ## Create vectors of the subject indices of the treatment units ordered by
  ## their matched set and corresponding control unit
  treated.subject.index=rep(0,max(matchedset.index.numeric.full,na.rm=T))
  matched.control.subject.index.mat=matrix(rep(0,nocontrols.per.match*length(treated.subject.index)),ncol=nocontrols.per.match)
  for(i in 1:length(treated.subject.index)){
    matched.set.temp=which(matchedset.index.numeric==i)
    treated.temp.index=which(datatemp$treated[matched.set.temp]==1)
    treated.subject.index[i]=matched.set.temp[treated.temp.index]
    matched.control.subject.index.mat[i,]=matched.set.temp[-treated.temp.index]
  }
  matched.control.subject.index=matched.control.subject.index.mat
  num_matched_sets=nrow(matched.control.subject.index.mat)
  ## Check balance
  missing.mat=matrix(rep(0,ncol(Xmat.full)*nrow(Xmat.full)),ncol=ncol(Xmat.full)
  )
  Xmat.without.missing=Xmat.full
  for(i in 1:ncol(Xmat.full)){
    Xmat.without.missing[missing.mat[,i]==1,i]=NA
  }
  ## Also compute balance on logit propensity score
  Xmat.without.missing$logit.ps=logit.propscore
  Xmat.without.missing=sapply(Xmat.without.missing,as.numeric)
  treatedmat=Xmat.without.missing[datatemp.full$treated==1,];
  ## Standardized differences before matching
  controlmat.before=Xmat.without.missing[datatemp.full$treated==0,];
  controlmean.before=apply(controlmat.before,2,mean,na.rm=TRUE);
  treatmean=apply(treatedmat,2,mean,na.rm=TRUE);
  treatvar=apply(treatedmat,2,var,na.rm=TRUE);
  controlvar=apply(controlmat.before,2,var,na.rm=TRUE);
  stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);
  ## Standardized differences after matching
  treatmat.after=Xmat.without.missing[treated.subject.index,]
  controlmat.after=Xmat.without.missing[matched.control.subject.index,];
  controlmean.after=apply(controlmat.after,2,mean,na.rm=TRUE);
  treatmean.after=apply(treatmat.after,2,mean,na.rm=TRUE)
  ## Standardized differences after matching
  stand.diff.after=(treatmean.after-controlmean.after)/sqrt((treatvar+controlvar)/2);
  sd.bf=stand.diff.before
  sd.af=stand.diff.after
  standardized_difs=cbind(sd.bf,sd.af)
  if(no.control.lack.overlap+no.treated.lack.overlap>0){
    data_outcomes <- rbind(data_outcomes[-which.remove,],data_outcomes[which.remove,])
  }
  ix_subjects_treated=treated.subject.index
  ix_subjects_control=matched.control.subject.index
  return(list(standardized_difs=standardized_difs,
              num_matched_sets=num_matched_sets,
              data_outcomes=data_outcomes,
              ix_subjects_treated=ix_subjects_treated,
              ix_subjects_control=ix_subjects_control))
}


filenames <- list.files(pattern='.XPT')
filelist <- vector("list",length=length(filenames))

for (fileidx in 1:length(filenames)){
  file <- read_xpt(filenames[fileidx])
  filelist[[fileidx]] <- file
}

dat <- filelist %>% reduce(full_join, by = "SEQN")

columnlabels <- map_df(dat, ~attributes(.x))$label

# treatment is "Ratio of family income to poverty" (col 147)
trtcolumn <- 147
columnlabels[trtcolumn]
dat <- dat[!is.na(dat[trtcolumn]),]
incometopov <- dat[trtcolumn]

treatment <- incometopov < 1
control <- incometopov > 1

trtvec <- rep(NA,length(incometopov))
trtvec[control] <- 0
trtvec[treatment] <- 1

dat <- dat[!is.na(trtvec),]
trtvec <- trtvec[!is.na(trtvec)]

trtvec <- trtvec[which(dat$RIDSTATR==2)] # must be interviewed and examined
dat <- dat[which(dat$RIDSTATR==2),] # must be interviewed and examined


# match on gender, race, age, etc.
matchcolumns <- c(133:136,138)
columnlabels[matchcolumns]
covariates <- dat[,matchcolumns]
# note factor variables
covariates[,c(1,3:5)] <- lapply(covariates[,c(1,3:5)],as.factor)


# matchres <- run_PSM(data_matching=covariates,data_outcomes=dat,treated=trtvec,calipersd=.5)
# matchres$standardized_difs
# plot.dataframe=data.frame(SMD=c(abs(matchres$standardized_difs[,2]),
#                                 abs(matchres$standardized_difs[,1])),
#                           Covariates=rep(rownames(matchres$standardized_difs),2),
#                           type=c(rep("After Matching",length(rownames(matchres$standardized_difs))),
#                                  rep("Before Matching",length(rownames(matchres$standardized_difs)))))
# ggplot(plot.dataframe,aes(x=SMD,y=Covariates))+
#   geom_point(size=3,aes(shape=factor(type,levels = c('Before Matching','After Matching'))))+
#   scale_shape_manual(values =c(21,16))+
#   scale_y_discrete(limits = rownames(matchres$standardized_difs))+
#   geom_vline(xintercept=c(0,0.1),lty=2) +
#   labs(x = "Absolute Standardized Mean Differences", y="Covariates") + theme_bw() +
#   theme(axis.text.y=element_text(size=8),legend.title = element_blank(),legend.position="bottom")
# 
# 
# 
# dfcov=df[,1:(ncol(df)-1)]
# calvec=rep(Inf,ncol(dfcov))
# names(calvec)<-names(dfcov)
# calvec[2]<-0.75
# m.out=MatchIt::matchit(trtvec~.,data=df,distance="robust_mahalanobis",caliper=calvec,std.caliper=T)
# summary(m.out)
# cobalt::bal.tab(m.out,m.threshold=.1)
# cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)


df=cbind(covariates,trtvec)
m.out=MatchIt::matchit(trtvec~.,data=df,method="optimal",distance="robust_mahalanobis")
summary(m.out)
cobalt::bal.tab(m.out,m.threshold=.1)
cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)


treated_idx=as.numeric(rownames(m.out$match.matrix))
control_idx=as.numeric(m.out$match.matrix)
num_matched_sets=length(treated_idx)



dat_treated <- dat[treated_idx,]
dat_control <- dat[control_idx,]


# divide data into planning and analysis samples
planning_sample_prop = 0.2
planning_sample_size_sets <- floor(planning_sample_prop*num_matched_sets)
planning_sample_size_subjects <- planning_sample_size_sets*2
ix_subjects_treated_planning <- sample(x=treated_idx,size=planning_sample_size_sets)
ix_subjects_control_planning <- control_idx[match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_treated_analysis <- treated_idx[-match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_control_analysis <- control_idx[match(ix_subjects_treated_analysis,treated_idx)]
# some relevant parameters
n1 <- planning_sample_size_sets
n2 <- num_matched_sets-planning_sample_size_sets


### EXAM
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(25,26,27,28,37,41,45,49,53,57,61,65,69,73,77,81,82,83,84,85,86,87,88)] # Cardioresp endurance
columnlabels[c(93,94,106,107,111,112,116,117,119,120)] # Cardiovascular Fitness
columnlabels[c(379:384)] # Lower body muscle strength
#columnlabels[c(seq(391,453,2),457,460:493)] # sports today or yesterday (technically part of lower body dataset)
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test
columnlabels[c(308:360)] # gross motor

## roughly 110


### QUESTIONNAIRE
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525:541,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition

## roughly 80

































































































































































































































































































































































































































































































































































































































































































































dat <- filelist %>% reduce(full_join, by = "SEQN")
columnlabels <- map_df(dat, ~attributes(.x))$label
# treatment is "Ratio of family income to poverty" (col 147)
trtcolumn <- 147
columnlabels[trtcolumn]
dat <- dat[!is.na(dat[trtcolumn]),]

dat <- dat[dat$RIDAGEYR >=3 & dat$RIDAGEYR <=5,]

incometopov <- dat[trtcolumn]
treatment <- incometopov < 1
control <- incometopov > 1
trtvec <- rep(NA,length(incometopov))
trtvec[control] <- 0
trtvec[treatment] <- 1
dat <- dat[!is.na(trtvec),]
trtvec <- trtvec[!is.na(trtvec)]
trtvec <- trtvec[which(dat$RIDSTATR==2)] # must be interviewed and examined
dat <- dat[which(dat$RIDSTATR==2),] # must be interviewed and examined
# match on gender, race, age, etc.
matchcolumns <- c(133:136,138)
columnlabels[matchcolumns]
covariates <- dat[,matchcolumns]
# note factor variables
covariates[,c(1,3:5)] <- lapply(covariates[,c(1,3:5)],as.factor)
df=cbind(covariates,trtvec)
#m.out=MatchIt::matchit(trtvec~.,data=df,distance = "glm", link = "probit",caliper = 2/3)#, exact = ~RIDRETH1)
m.out=MatchIt::matchit(trtvec~.,data=df,distance = "robust_mahalanobis", exact = ~RIDRETH1)
summary(m.out)
cobalt::bal.tab(m.out,m.threshold=.1)
cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)
treated_idx=as.numeric(rownames(m.out$match.matrix))
control_idx=as.numeric(m.out$match.matrix)
# Filter out pairs where either element is NA
filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
# Create new vectors without NA pairs
treated_idx <- treated_idx[filtered_indices]
control_idx <- control_idx[filtered_indices]
num_matched_sets=length(treated_idx)
dat_treated <- dat[treated_idx,]
dat_control <- dat[control_idx,]


columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(308:360)] # gross motor

## 61

columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition
columnlabels[c(525,544:577,580,581,587:609,613)] # phys activity
columnlabels[627] # general health condition

## realistically, less than 10



# divide data into planning and analysis samples
planning_sample_prop = 0.2
planning_sample_size_sets <- floor(planning_sample_prop*num_matched_sets)
planning_sample_size_subjects <- planning_sample_size_sets*2
ix_subjects_treated_planning <- sample(x=treated_idx,size=planning_sample_size_sets)
ix_subjects_control_planning <- control_idx[match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_treated_analysis <- treated_idx[-match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_control_analysis <- control_idx[match(ix_subjects_treated_analysis,treated_idx)]
# some relevant parameters
n1 <- planning_sample_size_sets
n2 <- num_matched_sets-planning_sample_size_sets



boxplot( dat[ix_subjects_treated_planning,c(3,7,11,13,15,17,19,21)]-dat[ix_subjects_control_planning,c(3,7,11,13,15,17,19,21)])
abline(h=0,col='red',lty=5)


boxplot( dat[ix_subjects_treated_planning,c(350:360)]-dat[ix_subjects_control_planning,c(350:360)])
abline(h=0,col='red',lty=5)


boxplot( dat[ix_subjects_treated_planning,c(302,303,373,627)]-dat[ix_subjects_control_planning,c(302,303,373,627)])
abline(h=0,col='red',lty=5)

boxplot( dat[ix_subjects_treated_planning,c(525,544:577,580,581,587:609,613)]-dat[ix_subjects_control_planning,c(525,544:577,580,581,587:609,613)])
abline(h=0,col='red',lty=5)


###
























dat <- filelist %>% reduce(full_join, by = "SEQN")
columnlabels <- map_df(dat, ~attributes(.x))$label
# treatment is "Ratio of family income to poverty" (col 147)
trtcolumn <- 147
columnlabels[trtcolumn]
dat <- dat[!is.na(dat[trtcolumn]),]

dat <- dat[dat$RIDAGEYR >=6 & dat$RIDAGEYR <=11,]

incometopov <- dat[trtcolumn]
treatment <- incometopov < 1
control <- incometopov > 1
trtvec <- rep(NA,length(incometopov))
trtvec[control] <- 0
trtvec[treatment] <- 1
dat <- dat[!is.na(trtvec),]
trtvec <- trtvec[!is.na(trtvec)]
trtvec <- trtvec[which(dat$RIDSTATR==2)] # must be interviewed and examined
dat <- dat[which(dat$RIDSTATR==2),] # must be interviewed and examined
# match on gender, race, age, etc.
matchcolumns <- c(133:136,138)
columnlabels[matchcolumns]
covariates <- dat[,matchcolumns]
# note factor variables
covariates[,c(1,3:5)] <- lapply(covariates[,c(1,3:5)],as.factor)
df=cbind(covariates,trtvec)
#m.out=MatchIt::matchit(trtvec~.,data=df,distance = "glm", link = "probit",caliper = 2/3)#, exact = ~RIDRETH1)
m.out=MatchIt::matchit(trtvec~.,data=df,distance = "robust_mahalanobis", exact = ~RIDRETH1)
summary(m.out)
cobalt::bal.tab(m.out,m.threshold=.1)
cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)
treated_idx=as.numeric(rownames(m.out$match.matrix))
control_idx=as.numeric(m.out$match.matrix)
# Filter out pairs where either element is NA
filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
# Create new vectors without NA pairs
treated_idx <- treated_idx[filtered_indices]
control_idx <- control_idx[filtered_indices]
num_matched_sets=length(treated_idx)
dat_treated <- dat[treated_idx,]
dat_control <- dat[control_idx,]



### EXAM
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(25,26,27,28,37,41,45,49,53,57,61,65,69,73,77,81,82,83,84,85,86,87,88)] # Cardioresp endurance
columnlabels[c(379:384)] # Lower body muscle strength
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test

## roughly 45


### QUESTIONNAIRE
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition

## realistically, less than 10























dat <- filelist %>% reduce(full_join, by = "SEQN")
columnlabels <- map_df(dat, ~attributes(.x))$label
# treatment is "Ratio of family income to poverty" (col 147)
trtcolumn <- 147
columnlabels[trtcolumn]
dat <- dat[!is.na(dat[trtcolumn]),]

dat <- dat[dat$RIDAGEYR >=12 & dat$RIDAGEYR <=15,]

incometopov <- dat[trtcolumn]
treatment <- incometopov < 1
control <- incometopov > 1
trtvec <- rep(NA,length(incometopov))
trtvec[control] <- 0
trtvec[treatment] <- 1
dat <- dat[!is.na(trtvec),]
trtvec <- trtvec[!is.na(trtvec)]
trtvec <- trtvec[which(dat$RIDSTATR==2)] # must be interviewed and examined
dat <- dat[which(dat$RIDSTATR==2),] # must be interviewed and examined
# match on gender, race, age, etc.
matchcolumns <- c(133:136,138)
columnlabels[matchcolumns]
covariates <- dat[,matchcolumns]
# note factor variables
covariates[,c(1,3:5)] <- lapply(covariates[,c(1,3:5)],as.factor)
df=cbind(covariates,trtvec)
#m.out=MatchIt::matchit(trtvec~.,data=df,distance = "glm", link = "probit",caliper = 2/3)#, exact = ~RIDRETH1)
m.out=MatchIt::matchit(trtvec~.,data=df,distance = "robust_mahalanobis", exact = ~RIDRETH1)
summary(m.out)
cobalt::bal.tab(m.out,m.threshold=.1)
cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)
treated_idx=as.numeric(rownames(m.out$match.matrix))
control_idx=as.numeric(m.out$match.matrix)
# Filter out pairs where either element is NA
filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
# Create new vectors without NA pairs
treated_idx <- treated_idx[filtered_indices]
control_idx <- control_idx[filtered_indices]
num_matched_sets=length(treated_idx)
dat_treated <- dat[treated_idx,]
dat_control <- dat[control_idx,]



### EXAM
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(93,94,106,107,111,112,116,117,119,120)] # Cardiovascular Fitness
columnlabels[c(379:384)] # Lower body muscle strength
#columnlabels[c(seq(391,453,2),457,460:493)] # sports today or yesterday (technically part of lower body dataset)
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test

## roughly 30


### QUESTIONNAIRE
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525:541,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition

## realistically 15ish


# divide data into planning and analysis samples
planning_sample_prop = 0.2
planning_sample_size_sets <- floor(planning_sample_prop*num_matched_sets)
planning_sample_size_subjects <- planning_sample_size_sets*2
ix_subjects_treated_planning <- sample(x=treated_idx,size=planning_sample_size_sets)
ix_subjects_control_planning <- control_idx[match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_treated_analysis <- treated_idx[-match(ix_subjects_treated_planning,treated_idx)]
ix_subjects_control_analysis <- control_idx[match(ix_subjects_treated_analysis,treated_idx)]
# some relevant parameters
n1 <- planning_sample_size_sets
n2 <- num_matched_sets-planning_sample_size_sets



datplot=dat[ix_subjects_treated_planning,c(7,93,94,521,524,626)]-dat[ix_subjects_control_planning,c(7,93,94,521,524,626)]
s=apply(datplot,2,function(x)sd(x,na.rm = T))
n=apply(datplot,2,function(x)sum(!is.na(x)))
bmiplot=datplot[,1]/s[1]
vo2plot=datplot[,2]/s[2]
cvfitplot=datplot[,3]/s[3]
gripplot=datplot[,4]/s[4]
pullplot=datplot[,5]/s[5]
plankplot=datplot[,6]/s[6]

boxplot( list(bmiplot,vo2plot,cvfitplot,gripplot,pullplot,plankplot) )
abline(h=0,col='red',lty=5)


dfplot=as.data.frame(list(bmiplot,vo2plot,cvfitplot,gripplot,pullplot,plankplot))
names(dfplot)=c('BMI','VO2MAX', 'ENDURANCE', 'GRIP', 'PULLUP', 'PLANK')

ggplot(reshape2::melt(dfplot), aes(x=variable, y=value)) +
  geom_boxplot(fill=rgb(232/255, 244/255, 255/255)) + theme_bw() + labs(x='',y='Poor - Not Poor')+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)+ylim(-3,3)+ theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())

