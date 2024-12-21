
##### load and process data #####
setwd("C:/Users/bekerman/Downloads/relevant-xpt-files")
library(tidyverse)
library(MASS)
library(haven)
library(caret)
library(optmatch)

filenames <- list.files(pattern='.XPT')
filelist <- vector("list",length=length(filenames))
for (fileidx in 1:length(filenames)){
  file <- read_xpt(filenames[fileidx])
  filelist[[fileidx]] <- file
}
dat <- filelist %>% reduce(full_join, by = "SEQN")
dat <- dat %>% group_by(SEQN) %>% slice(1) # resolve issue w/ multiples of same row
columnlabels <- map_df(dat, ~attributes(.x))$label

### EXAM (~ 100 outcomes)
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(25,26,27,28,37,41,45,49,53,57,61,65,69,73,77,81,82,83,84,85,86,87,88)] # Cardioresp endurance
columnlabels[c(93,94,106,107,111,112,116,117,119,120)] # Cardiovascular Fitness
columnlabels[c(379:384)] # Lower body muscle strength
#columnlabels[c(seq(391,453,2),457,460:493)] # sports today or yesterday (technically part of lower body dataset)
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test
columnlabels[c(308:360)] # gross motor

### QUESTIONNAIRE (~ 80 outcomes)
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525:541,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition



##### define functions to get + plot data for age groups #####
# treatment is "ratio of family income to poverty" (FIPR; column 147)
trtcolumn <- 147
columnlabels[trtcolumn]
# match on gender, race, age, etc.
matchcolumns <- c(133:136,138)
columnlabels[matchcolumns]
# set planning/analysis ratio
planning_sample_prop = 0.2

# from exploration, change some columns in dat
zeros_PAQ605 <- dat$PAQ605==2
dat$PAQ610[zeros_PAQ605] <- 0; dat$PAD615[zeros_PAQ605] <- 0
zeros_PAQ620 <- dat$PAQ620==2
dat$PAQ625[zeros_PAQ620] <- 0; dat$PAQ625[dat$PAQ625==99] <- NA; dat$PAD630[zeros_PAQ620] <- 0
zeros_PAQ635 <- dat$PAQ635==2
dat$PAQ640[zeros_PAQ635] <- 0; dat$PAD645[zeros_PAQ635] <- 0
zeros_PAQ650 <- dat$PAQ650==2
dat$PAQ655[zeros_PAQ650] <- 0; dat$PAD660[zeros_PAQ650] <- 0
zeros_PAQ665 <- dat$PAQ665==2
dat$PAQ670[zeros_PAQ665] <- 0; dat$PAD675[zeros_PAQ665] <- 0

dat$PHYSACTSUM <- apply( dat[,545:577],1,function(x)sum(x>0&x!=77&x!=99,na.rm=T) )
dat$PHYSACTSUM[is.na(dat$PAQ722)|dat$PAQ722==9] <- NA
dat$PARTICPSUM <- apply( dat[,587:609],1,function(x)sum(x>0&x!=77&x!=99,na.rm=T) )
dat$PARTICPSUM[is.na(dat$PAQ755)|dat$PAQ755==9] <- NA
dat$CARDIO_END_PROBS <- apply( dat[,83:88],1,function(x)sum(x>0,na.rm=T) )
dat$CARDIO_END_PROBS[dat$CEDEXSTS==2|dat$CEDEXSTS==3] <- NA

dat <- dat %>% mutate(minperwk_vigorouswork=prod(PAQ610,PAD615,na.rm=T),
                      minperwk_moderatework=prod(PAQ625,PAD630,na.rm=T),
                      minperwk_walkbike=prod(PAQ640,PAD645,na.rm=T),
                      minperwk_vigorousrec=prod(PAQ655,PAD660,na.rm=T),
                      minperwk_moderaterec=prod(PAQ670,PAD675,na.rm=T))
dat$minperwk_vigorouswork[is.na(dat$PAQ610)|is.na(dat$PAD615)] <- NA
dat$minperwk_moderatework[is.na(dat$PAQ625)|is.na(dat$PAD630)] <- NA
dat$minperwk_walkbike[is.na(dat$PAQ640)|is.na(dat$PAD645)] <- NA
dat$minperwk_vigorousrec[is.na(dat$PAQ655)|is.na(dat$PAD660)] <- NA
dat$minperwk_moderaterec[is.na(dat$PAQ670)|is.na(dat$PAD675)] <- NA

columnlabels <- c(columnlabels,"PHYSACTSUM","PARTICPSUM","CARDIO_END_PROBS",
                  "minperwk_vigorouswork","minperwk_moderatework",
                  "minperwk_walkbike","minperwk_vigorousrec","minperwk_moderaterec")

# define functions
get_age_data <- function(age_lower, age_upper, fipr_cutoff=1, files=filelist,
                         onlyPS=F,link,cal,useglmnet=F,nfolds=NULL,
                         num_controls=1,exact=NULL){
  # treatment is "Ratio of family income to poverty" (col 147)
  datold <- dat
  dat <- dat[!is.na(dat[trtcolumn]),]
  dat <- dat[dat$RIDAGEYR >=age_lower & dat$RIDAGEYR <=age_upper,]
  incometopov <- dat[trtcolumn]
  treatment <- incometopov < fipr_cutoff
  control <- incometopov > fipr_cutoff
  trtvec <- rep(NA,length(incometopov))
  trtvec[control] <- 0
  trtvec[treatment] <- 1
  dat <- dat[!is.na(trtvec),]
  trtvec <- trtvec[!is.na(trtvec)]
  trtvec <- trtvec[which(dat$RIDSTATR==2)] # must be interviewed and examined
  dat <- dat[which(dat$RIDSTATR==2),] # must be interviewed and examined
  # match on gender, race, age, etc.
  covariates <- dat[,matchcolumns]
  # note factor variables
  covariates[,c(1,3:5)] <- lapply(covariates[,c(1,3:5)],as.factor)
  covariates <- covariates %>% mutate(BLACK = ifelse(RIDRETH1==4, 1, 0))
  covariates <- covariates %>% mutate(WHITE = ifelse(RIDRETH1==3, 1, 0))
  # covariates <- covariates %>% mutate(HISP = ifelse(RIDRETH1 %in% c(1,2), 1, 0))
  covariates <- covariates %>% mutate(OTHER = ifelse(RIDRETH1 %in% c(1,2,5), 1, 0))
  covariates[,c(6:8)] <- lapply(covariates[,c(6:8)],as.factor)
  covariates <- subset(covariates,select=-c(RIDRETH1))
  colnames(covariates) <- c('Gender', 'Age', 'Time Period', 'Birth Country',
                            'Black', 'White', 'Other')
  df=cbind(covariates,trtvec)
  ##m.out=MatchIt::matchit(trtvec~.,data=df,distance = "robust_mahalanobis", exact = ~RIDRETH1)
  # m.out=MatchIt::matchit(trtvec~.,data=df,distance = "glm", link = "probit",caliper = 3/4)
  # summary(m.out)
  # bal_tbl <- cobalt::bal.tab(m.out,m.threshold=.1)
  # love_plt <- cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)
  # treated_idx=as.numeric(rownames(m.out$match.matrix))
  # control_idx=as.numeric(m.out$match.matrix)
  # # Filter out pairs where either element is NA
  # filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
  # # Create new vectors without NA pairs
  # treated_idx <- treated_idx[filtered_indices]
  # control_idx <- control_idx[filtered_indices]
  # num_matched_sets <- length(treated_idx)
  # return(list(matchvec=res$matchvec,
  #             which.remove=res$which.remove,
  #             treated_idx=res$treated_idx,
  #             control_idx=res$control_idx,
  #             num_matched_sets=length(res$treated_idx),
  #             standardized_difs=res$standardized_difs,
  #             love_plt=res$plt))
  
  data_matching=covariates
  treated=trtvec
  dff=cbind(data_matching,treated)
  
  ## internal matching function
  fxx<-function(onlyPS,link,cal,useglmnet,nfolds,num_controls,exact){
    # browser()
    ## Load libraries and auxillary functions
    library(MASS)
    library(sensitivitymult)
    library(optmatch)
    options("optmatch_max_problem_size" = Inf)
    ## Function for computing
    ## rank based Mahalanobis distance. Prevents an outlier from
    ## inflating the variance for a variable, thereby decreasing its importance.
    ## Also, the variances are not permitted to decrease as ties
    ## become more common, so that, for example, it is not more important
    ## to match on a rare binary variable than on a common binary variable
    ## z is a vector, length(z)=n, with z=1 for treated, z=0 for control
    ## X is a matrix with n rows containing variables in the distance
    ## onlYPS is boolean indicating whether to only match on propensity score dist.
    ## exact is, if non-null, the column indices of X on which to exactly match.
    smahal=function(z,X,onlyPS=F,exact=NULL,penalty=1000000){
      if (onlyPS & !is.null(exact)){
        Xwrong<-as.matrix(X)
        n<-dim(Xwrong)[1]
        rownames(Xwrong)<-1:n
        k<-dim(Xwrong)[2]
        m<-sum(z)
        for (j in 1:k) Xwrong[,j]<-rank(Xwrong[,j])
        cv<-cov(Xwrong)
        vuntied<-var(1:n)
        rat<-sqrt(vuntied/diag(cv))
        cv<-diag(rat)%*%cv%*%diag(rat)
        outwrong<-matrix(NA,m,n-m)
        Xc<-as.matrix(Xwrong[z==0,])
        Xt<-as.matrix(X[z==1,])
        rownames(outwrong)<-rownames(Xwrong)[z==1]
        colnames(outwrong)<-rownames(Xwrong)[z==0]
        icov<-ginv(cv)
        for (i in 1:m) outwrong[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
        for (i in 1:m){
          for (ex in exact){
            notequals <- Xt[i,ex]!=Xc[,ex]
            if (any(notequals)) outwrong[i,which(notequals)]=NA
          }
        }
        X=X$logit.ps
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
        out[which(is.na(outwrong))] <- penalty
      } else{
        if (onlyPS) X=X$logit.ps
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
        if (!is.null(exact)){
          for (i in 1:m){
            for (ex in exact){
              notequals <- Xt[i,ex]!=Xc[,ex]
              if (any(notequals)) outwrong[i,which(notequals)]=penalty
            }
          }
        }
      }
      out
    }
    ## Function for adding a propensity score caliper to a distance matrix dmat
    ## calipersd is the caliper in terms of standard deviation of the logit propensity score
    addcaliper=function(dmat,z,logitp,calipersd=.5,penalty=1000){
      # Pooled within group standard devation
      sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
      adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
      adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
      dmat=dmat+adif*penalty
      dmat
    }
    
    ## Propensity score model
    dff$treated=as.numeric(as.character(dff$treated))
    datatemp=dff
    if(useglmnet){
      library(glmnet)
      best_lambda = cv.glmnet(x=model.matrix( ~.,subset(datatemp, select=-c(treated)) )[,-1], y=datatemp$treated, family=binomial(link=link), nfolds=nfolds, alpha=0)$lambda.min
      propscore.model = glmnet(x=model.matrix( ~.,subset(datatemp, select=-c(treated)) )[,-1], y=datatemp$treated, family=binomial(link=link), alpha=0, lambda = best_lambda)
      dmy=dummyVars('treated~.',data=datatemp,fullRank=T)
      Xmat=data.frame(predict(dmy,newdata=datatemp))
      Xmatmahal=Xmat
      datatemp$logit.ps=as.numeric(predict(propscore.model,newx=as.matrix(Xmat)))
    } else{
      propscore.model = glm(treated~., family=binomial(link=link), data=datatemp)
      dmy=dummyVars('treated~.',data=datatemp)
      Xmat=data.frame(predict(dmy,newdata=datatemp))
      Xmatmahal=Xmat
      datatemp$logit.ps=as.numeric(predict(propscore.model,newx=as.matrix(Xmat)))
    }
    ## Use Hansen (2009)’s rule for removing subjects who lack overlap
    logit.propscore=datatemp$logit.ps
    pooled.sd.logit.propscore=sqrt(var(logit.propscore[datatemp$treated==1])/2+var
                                   (logit.propscore[datatemp$treated==0])/2)
    min.treated.logit.propscore=min(logit.propscore[datatemp$treated==1])
    max.control.logit.propscore=max(logit.propscore[datatemp$treated==0])
    ## How many treated and control subjects lack overlap by Hansen's criterion
    no.treated.lack.overlap=sum(logit.propscore[datatemp$treated==1]>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))
    no.control.lack.overlap=sum(logit.propscore[datatemp$treated==0]<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore))
    ## If there are subjects who lack overlap, remove them from the datatemp dataset
    datatemp.original=datatemp
    datatemp.full=datatemp
    Xmat.original=Xmat
    Xmat.full=Xmat
    which.remove=NA
    if(no.treated.lack.overlap+no.control.lack.overlap>0){
      which.remove=which((logit.propscore>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))|(logit.propscore<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore)))
      datatemp=datatemp[-which.remove,]
      datatemp.full=rbind(datatemp,datatemp.original[which.remove,])
      Xmat=Xmat[-which.remove,]
      Xmat.full=rbind(Xmat,Xmat.original[which.remove,])
      Xmatmahal=Xmatmahal[-which.remove,]
      logit.propscore=c(logit.propscore[-which.remove],logit.propscore[which.remove])
    }
    ## For the purposes of balance checking later, in datatemp.full, append
    ## the removed rows of datatemp to the end of datatemp
    rownames(datatemp)=seq(1,nrow(datatemp),1)
    Xmatmahal$logit.ps=datatemp$logit.ps
    ## Rank based Mahalanobis distance with caliper
    distmat=smahal(datatemp$treated,Xmatmahal,onlyPS=onlyPS,exact=exact)
    distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=cal)
    rownames(distmat)=rownames(datatemp)[datatemp$treated==1]
    colnames(distmat)=rownames(datatemp)[datatemp$treated==0]
    ## Match nocontrols.per.match to each treated unit
    nocontrols.per.match=num_controls
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
    treated_idx=treated.subject.index
    control_idx=matched.control.subject.index
    # # Filter out pairs where either element is NA
    # filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
    # # Create new vectors without NA pairs
    # treated_idx <- treated_idx[filtered_indices]
    # control_idx <- control_idx[filtered_indices]
    # num_matched_sets <- length(treated_idx)
    if (!is.na(which.remove)){
      data_matching <- rbind(data_matching[-which.remove,],data_matching[which.remove,])
    }
    ## get balance
    Xmat.without.missing=data_matching
    dmy=dummyVars(cbind(1,Xmat.without.missing),data=cbind(1,Xmat.without.missing),fullRank=T)
    Xmat.without.missing=data.frame(predict(dmy,newdata=cbind(1,Xmat.without.missing)))
    ## Also compute balance on logit propensity score
    #Xmat.without.missing$logit.ps=logit.propscore
    treatedmat=Xmat.without.missing[datatemp.full$treated==1,];
    ## Raw differences for binary variables
    binaryvars <- apply(Xmat.without.missing,2,function(x)all(x%in%(0:1)))
    ## Standardized differences before matching
    controlmat.before=Xmat.without.missing[datatemp.full$treated==0,];
    controlmean.before=apply(controlmat.before,2,mean,na.rm=TRUE);
    treatmean=apply(treatedmat,2,mean,na.rm=TRUE);
    treatvar=apply(treatedmat,2,var,na.rm=TRUE);
    controlvar=apply(controlmat.before,2,var,na.rm=TRUE);
    stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);
    stand.diff.before[binaryvars]=(treatmean-controlmean.before)[binaryvars]
    ## Standardized differences after matching
    treatmat.after=Xmat.without.missing[treated_idx,]
    controlmat.after=Xmat.without.missing[control_idx,];
    controlmean.after=apply(controlmat.after,2,mean,na.rm=TRUE);
    treatmean.after=apply(treatmat.after,2,mean,na.rm=TRUE)
    ## Standardized differences after matching
    stand.diff.after=(treatmean.after-controlmean.after)/sqrt((treatvar+controlvar)/2);
    stand.diff.after[binaryvars]=(treatmean.after-controlmean.after)[binaryvars]
    # sd.bf=stand.diff.before[1:(length(stand.diff.before)-1)]
    # sd.af=stand.diff.after[1:(length(stand.diff.after)-1)]
    sd.bf=stand.diff.before[1:(length(stand.diff.before))]
    sd.af=stand.diff.after[1:(length(stand.diff.after))]
    standardized_difs=cbind(sd.bf,sd.af)
    rownames(standardized_difs) <- c('Gender', 'Age', 'Time Period', 'Birth Country',
                                     'Black', 'White', 'Other')
    plot.dataframe=data.frame(SMD=c(abs(sd.af),abs(sd.bf)),Covariates=rep(names(sd.bf),2),type=c(rep("After Matching",length(names(sd.bf))),rep("Before Matching",length(names(sd.bf)))))
    cov_order=names(sd.af)
    plt=ggplot(plot.dataframe,aes(x=SMD,y=Covariates))+
      geom_point(size=3,aes(shape=factor(type,levels = c('Before Matching','After Matching'))))+
      scale_shape_manual(values =c(21,16))+
      scale_y_discrete(limits = rev(cov_order), labels=rev(c('Gender', 'Age', 'Time Period', 'Birth Country',
                                                    'Black', 'White', 'Other')))+
      geom_vline(xintercept=c(0,0.1),lty=2) +
      xlim(0,.4) +
      labs(x = "Absolute Differences", y="Covariates") + theme_bw() +
      theme(axis.text.y=element_text(size=11),
            axis.text.x=element_text(size=11),
            axis.title.x=element_text(size=11),
            axis.title.y=element_text(size=11),
            legend.text=element_text(size=11),
            legend.title = element_blank(),legend.position="bottom")
    return(list(matchvec=matchvec,
                which.remove=which.remove,
                treated_idx=treated_idx,
                control_idx=control_idx,
                num_matched_sets=length(treated_idx),
                standardized_difs=standardized_difs,
                plt=plt))
  }
  res=fxx(onlyPS=onlyPS,link=link,cal=cal,useglmnet=useglmnet,nfolds=nfolds,num_controls=num_controls,exact=exact)
  datnew <- dat; dat <- datold
  return(list(dat=datnew,
              treated_idx=res$treated_idx,
              control_idx=res$control_idx,
              num_matched_sets=length(res$treated_idx),
              bal_tbl=res$standardized_difs,
              love_plt=res$plt))
}
plot_age_data <- function(age_data, treated, control, col_num){
  tminusc <- age_data[treated,col_num]-age_data[control,col_num]
  s <- apply(as.matrix(tminusc),2,function(x)sd(x,na.rm=T))
  n <- apply(as.matrix(tminusc),2,function(x)sum(!is.na(x)))
  dat_plt <- tminusc#/(s/sqrt(n))
  
  dfplot=as.data.frame( list(dat_plt) )
  names(dfplot)=columnlabels[col_num]
  
  # plt <- ggplot(reshape2::melt(dfplot), aes(x=variable, y=value)) +
  #   geom_boxplot(fill=rgb(232/255, 244/255, 255/255)) + theme_bw() + ggtitle(paste('t-Scores:',columnlabels[col_num])) + labs(x='',y='Poor - Not Poor')+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)+ theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())#+ylim(-3,3)+ theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
  
  plt <- ggplot(reshape2::melt(dfplot), aes(x=variable, y=value)) +
    geom_boxplot(fill=rgb(232/255, 244/255, 255/255)) + theme_bw() + ggtitle(paste(columnlabels[col_num])) + labs(x='',y='Poor - Not Poor')+geom_hline(yintercept=0,linetype='dotted',color=rgb(39/255, 60/255, 117/255),linewidth=1.3)+ theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())#+ylim(-3,3)+ theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
  
  return(plt)
}


##### Ages 3-5 #####
set.seed(0) # for reproducibility
res_3to5 <- get_age_data(age_lower=3, age_upper=5, link='probit', cal=0.5)
dat_3to5 <- res_3to5$dat
treated_idx_3to5 <- res_3to5$treated_idx
control_idx_3to5 <- res_3to5$control_idx
num_matched_sets_3to5 <- res_3to5$num_matched_sets

### EXAM (~ 60 outcomes)
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(308:360)] # gross motor

### QUESTIONNAIRE (< 10 outcomes)
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition
columnlabels[c(525,544:577,580,581,587:609,613)] # phys activity
columnlabels[627] # general health condition

# divide data into planning and analysis samples
all_colnums_3to5 <- c(3,7,9,11,13,15,17,19,21,308:355,358,359,303,373,525,
                      627,638,639)
planning_sample_size_sets_3to5 <- floor(planning_sample_prop*num_matched_sets_3to5)
planning_sample_size_subjects_3to5 <- planning_sample_size_sets_3to5*2
ix_subjects_treated_planning_3to5 <- sample(x=treated_idx_3to5,size=planning_sample_size_sets_3to5)
ix_subjects_control_planning_3to5 <- control_idx_3to5[match(ix_subjects_treated_planning_3to5,treated_idx_3to5)]
ix_subjects_treated_analysis_3to5 <- treated_idx_3to5[-match(ix_subjects_treated_planning_3to5,treated_idx_3to5)]
ix_subjects_control_analysis_3to5 <- control_idx_3to5[match(ix_subjects_treated_analysis_3to5,treated_idx_3to5)]
n1 <- planning_sample_size_sets_3to5
n2 <- num_matched_sets_3to5-n1



##### Ages 6-11 #####
set.seed(0) # for reproducibility
res_6to11 <- get_age_data(age_lower=6, age_upper=11, link='probit', cal=0.5)
dat_6to11 <- res_6to11$dat
treated_idx_6to11 <- res_6to11$treated_idx
control_idx_6to11 <- res_6to11$control_idx
num_matched_sets_6to11 <- res_6to11$num_matched_sets

### EXAM (~ 45 outcomes)
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(25,26,27,28,37,41,45,49,53,57,61,65,69,73,77,81,82,83,84,85,86,87,88)] # Cardioresp endurance
columnlabels[c(379:384)] # Lower body muscle strength
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test

### QUESTIONNAIRE (< 10 outcomes)
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition

# divide data into planning and analysis samples
all_colnums_6to11 <- c(3,7,9,11,13,15,17,19,21,25,26,27,28,37,41,45,49,53,57,61,65,
                       69,73,77,81,82,379:384,524,509,511,513,
                       515,517,519,303,373,525,
                       615,627,638,639,640)
planning_sample_size_sets_6to11 <- floor(planning_sample_prop*num_matched_sets_6to11)
planning_sample_size_subjects_6to11 <- planning_sample_size_sets_6to11*2
ix_subjects_treated_planning_6to11 <- sample(x=treated_idx_6to11,size=planning_sample_size_sets_6to11)
ix_subjects_control_planning_6to11 <- control_idx_6to11[match(ix_subjects_treated_planning_6to11,treated_idx_6to11)]
ix_subjects_treated_analysis_6to11 <- treated_idx_6to11[-match(ix_subjects_treated_planning_6to11,treated_idx_6to11)]
ix_subjects_control_analysis_6to11 <- control_idx_6to11[match(ix_subjects_treated_analysis_6to11,treated_idx_6to11)]
n1 <- planning_sample_size_sets_6to11
n2 <- num_matched_sets_6to11-n1



##### Ages 12-15 #####
set.seed(0) # for reproducibility
res_12to15 <- get_age_data(age_lower=12, age_upper=15, link='probit', cal=0.5)
dat_12to15 <- res_12to15$dat
treated_idx_12to15 <- res_12to15$treated_idx
control_idx_12to15 <- res_12to15$control_idx
num_matched_sets_12to15 <- res_12to15$num_matched_sets

### EXAM (~ 30 outcomes)
columnlabels[c(3,7,11,13,15,17,19,21)] # body measures
columnlabels[c(93,94,106,107,111,112,116,117,119,120)] # Cardiovascular Fitness
columnlabels[c(379:384)] # Lower body muscle strength
#columnlabels[c(seq(391,453,2),457,460:493)] # sports today or yesterday (technically part of lower body dataset)
columnlabels[c(524)] # modified pullup
columnlabels[c(509,511,513,515,517,519,521)] # grip test

### QUESTIONNAIRE (~ 15 outcomes)
columnlabels[c(302,303)] # child overweight
columnlabels[373] # general health condition    
columnlabels[c(525:541,544:577,580,581,587:609,613,615)] # phys activity
columnlabels[627] # general health condition

# divide data into planning and analysis samples
all_colnums_12to15 <- c(3,7,9,11,13,15,17,19,21,93,94,106,107,111,112,116,117,119,120,379:384,
                        524,509,511,513,515,517,519,303,373,525,580,581,
                        613,615,627,638,639,641:645)
planning_sample_size_sets_12to15 <- floor(planning_sample_prop*num_matched_sets_12to15)
planning_sample_size_subjects_12to15 <- planning_sample_size_sets_12to15*2
ix_subjects_treated_planning_12to15 <- sample(x=treated_idx_12to15,size=planning_sample_size_sets_12to15)
ix_subjects_control_planning_12to15 <- control_idx_12to15[match(ix_subjects_treated_planning_12to15,treated_idx_12to15)]
ix_subjects_treated_analysis_12to15 <- treated_idx_12to15[-match(ix_subjects_treated_planning_12to15,treated_idx_12to15)]
ix_subjects_control_analysis_12to15 <- control_idx_12to15[match(ix_subjects_treated_analysis_12to15,treated_idx_12to15)]
n1 <- planning_sample_size_sets_12to15
n2 <- num_matched_sets_12to15-n1



##### exploration on planning sample #####
tbl_3to5 <- apply(dat_3to5[,all_colnums_3to5],2,function(x) table(x,useNA = 'always')) #65 outcomes
names(tbl_3to5) <- columnlabels[all_colnums_3to5]
tbl_3to5
tbl_6to11 <- apply(dat_6to11[,all_colnums_6to11],2,function(x) table(x,useNA = 'always')) #47 outcomes
names(tbl_6to11) <- columnlabels[all_colnums_6to11]
tbl_6to11
tbl_12to15 <- apply(dat_12to15[,all_colnums_12to15],2,function(x) table(x,useNA = 'always')) #47 outcomes
names(tbl_12to15) <- columnlabels[all_colnums_12to15]
tbl_12to15

# plot_age_data(dat_12to15, ix_subjects_treated_planning_12to15, ix_subjects_control_planning_12to15, col_num=7)

res_3to5$love_plt; res_3to5$num_matched_sets
res_6to11$love_plt; res_6to11$num_matched_sets
res_12to15$love_plt; res_12to15$num_matched_sets
