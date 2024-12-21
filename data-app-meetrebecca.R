
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
columnlabels <- c(columnlabels,"PHYSACTSUM","PARTICPSUM","CARDIO_END_PROBS")

# define functions
get_age_data <- function(age_lower, age_upper, fipr_cutoff=1, files=filelist){
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
  df=cbind(covariates,trtvec)
  #m.out=MatchIt::matchit(trtvec~.,data=df,distance = "robust_mahalanobis", exact = ~RIDRETH1)
  m.out=MatchIt::matchit(trtvec~.,data=df,distance = "glm", link = "probit",caliper = 3/4)
  summary(m.out)
  bal_tbl <- cobalt::bal.tab(m.out,m.threshold=.1)
  love_plt <- cobalt::love.plot(cobalt::bal.tab(m.out,m.threshold=.1),stars = 'std', abs=T)
  treated_idx=as.numeric(rownames(m.out$match.matrix))
  control_idx=as.numeric(m.out$match.matrix)
  # Filter out pairs where either element is NA
  filtered_indices <- !is.na(treated_idx) & !is.na(control_idx)
  # Create new vectors without NA pairs
  treated_idx <- treated_idx[filtered_indices]
  control_idx <- control_idx[filtered_indices]
  num_matched_sets <- length(treated_idx)
  datnew <- dat; dat <- datold
  return(list(dat=datnew,
              treated_idx=treated_idx,
              control_idx=control_idx,
              num_matched_sets=num_matched_sets,
              bal_tbl=bal_tbl,
              love_plt=love_plt))
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
res_3to5 <- get_age_data(age_lower=3, age_upper=5)
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
all_colnums_3to5 <- c(3,7,9,11,13,15,17,19,21,308:355,358,359,302,303,373,525,
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
res_6to11 <- get_age_data(age_lower=6, age_upper=11)
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
                       515,517,519,302,303,373,525,
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
res_12to15 <- get_age_data(age_lower=12, age_upper=15)
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
                        524,509,511,513,515,517,519,302,303,373,(525:541)[-c(2,5,8,11,14)],580,581,
                        613,615,627,638,639)
planning_sample_size_sets_12to15 <- floor(planning_sample_prop*num_matched_sets_12to15)
planning_sample_size_subjects_12to15 <- planning_sample_size_sets_12to15*2
ix_subjects_treated_planning_12to15 <- sample(x=treated_idx_12to15,size=planning_sample_size_sets_12to15)
ix_subjects_control_planning_12to15 <- control_idx_12to15[match(ix_subjects_treated_planning_12to15,treated_idx_12to15)]
ix_subjects_treated_analysis_12to15 <- treated_idx_12to15[-match(ix_subjects_treated_planning_12to15,treated_idx_12to15)]
ix_subjects_control_analysis_12to15 <- control_idx_12to15[match(ix_subjects_treated_analysis_12to15,treated_idx_12to15)]
n1 <- planning_sample_size_sets_12to15
n2 <- num_matched_sets_12to15-n1



##### exploration on planning sample #####
tbl_3to5 <- apply(dat_3to5[,all_colnums_3to5],2,function(x) table(x,useNA = 'always')) #66 outcomes
names(tbl_3to5) <- columnlabels[all_colnums_3to5]
tbl_3to5
tbl_6to11 <- apply(dat_6to11[,all_colnums_6to11],2,function(x) table(x,useNA = 'always')) #48 outcomes
names(tbl_6to11) <- columnlabels[all_colnums_6to11]
tbl_6to11
tbl_12to15 <- apply(dat_12to15[,all_colnums_12to15],2,function(x) table(x,useNA = 'always')) #54 outcomes
names(tbl_12to15) <- columnlabels[all_colnums_12to15]
tbl_12to15

# plot_age_data(dat_12to15, ix_subjects_treated_planning_12to15, ix_subjects_control_planning_12to15, col_num=7)#94)
  
res_3to5$love_plt; res_3to5$num_matched_sets
res_6to11$love_plt; res_6to11$num_matched_sets
res_12to15$love_plt; res_12to15$num_matched_sets
