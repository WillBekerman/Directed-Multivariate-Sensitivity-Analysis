
################################################################################
#  Load and process data
################################################################################
library(tidyverse)
library(MASS)
library(haven)
library(caret)
library(optmatch)
library(gridExtra)
library(purrr)
library(MatchIt)
library(stringr)

dir <- "nhanes_data"
filenames <- list.files(dir, pattern='.xpt')
filelist <- vector("list",length=length(filenames))
for (fileidx in 1:length(filenames)){
  file <- read_xpt(paste0(dir,'/',filenames[fileidx]))
  filelist[[fileidx]] <- file
}

# Add wave record to each file
files_tbl <- tibble(file = filenames) %>%
  mutate(
    WAVE = case_when(
      TRUE ~ str_match(file, "_([A-Z])\\.xpt$")[,2]
    ),
    dataset = case_when(
      TRUE ~ str_remove(file, "_[A-Z]\\.xpt$")
    ),
    WAVE = factor(WAVE, levels = c("A","B","C","D","E",
                                   "F","G","H","I"))
)
waves_chr <- tibble(file = filenames) %>%
  mutate(WAVE = str_match(file, "_([A-Z])\\.xpt$")[,2]) %>%
  pull(WAVE)
wave_levels <- c("A","B","C","D","E","F","G","H","I")
filelist <- map2(
  filelist, waves_chr,
  ~ mutate(.x, WAVE = factor(.y, levels = wave_levels, labels = 1:9))
)

# Functions to combine files properly into same column names
# Safely read a label from common places
get_label <- function(x) {
  lab <- attr(x, "label", exact = TRUE)
  if (is.null(lab)) "" else as.character(lab)
}
# Reduce full_join, coalesce duplicate columns, and preserve labels
full_join_reduce_coalesce_clean <- function(df_list, by, suffix = c("", ".dup")) {
  stopifnot(length(df_list) >= 1)
  out <- reduce(df_list, ~ full_join(.x, .y, by = by, suffix = suffix))
  nms  <- names(out)
  base <- sub("(?:\\.(?:x|y)|\\.dup)+$", "", nms, perl = TRUE)
  groups <- split(seq_along(nms), base)
  for (bn in names(groups)) {
    idx  <- groups[[bn]]
    cols <- nms[idx]
    if (length(cols) > 1) {
      labels <- map_chr(cols, ~ get_label(out[[.x]]))
      winner <- labels[which(nzchar(labels))[1]] %||% ""
      out[[bn]] <- dplyr::coalesce(!!!out[cols])
      if (nzchar(winner)) attr(out[[bn]], "label") <- winner
      out[setdiff(cols, bn)] <- NULL
    }
  }
  names(out) <- sub("(?:\\.(?:x|y)|\\.dup)+$", "", names(out), perl = TRUE)
  out
}
column_labels <- function(dat, default = c("name", "empty")) {
  default <- match.arg(default)
  labs <- map_chr(dat, get_label)
  if (default == "name") {
    labs <- ifelse(nzchar(labs), labs, names(dat))
  } else {
    labs <- ifelse(nzchar(labs), labs, "")
  }
  setNames(labs, names(dat))
}

# Process early waves of physical activity separately
files_iaf <- unlist(lapply(filelist,function(x)
  'PADLEVEL'%in%colnames(x))) # get PAQIAF files
files_iaf_paq <- which(unlist(lapply(filelist,function(x)
  'PAD200'%in%colnames(x)))) # get related PAQ files
filelist[files_iaf_paq] <- lapply(filelist[files_iaf_paq], function(df) {
  df %>%
    mutate(
      PAD200 = case_when(
        PAD200 %in% c(2, 3) ~ 0,
        PAD200 %in% c(7, 7777, 9, 9999) ~ NA_real_,
        TRUE ~ PAD200
      ),
      PAD320 = case_when(
        PAD320 %in% c(2, 3) ~ 0,
        PAD320 %in% c(7, 7777, 9, 9999) ~ NA_real_,
        TRUE ~ PAD320
      )
    )
})
paq_df <- bind_rows(filelist[files_iaf_paq])  # has PAD200, PAD320
iaf_df <- bind_rows(filelist[files_iaf])      # has PADLEVEL, PADTIMES, PADDURAT
merged_df <- full_join(paq_df, iaf_df, by = "SEQN")
new_df <- merged_df %>%
  group_by(SEQN) %>%
  summarise(
    mins_vigorous_act_wk_earlywaves = case_when(
      all(PAD200 == 0) ~ 0,  # PAQ says no vigorous
      any(PADLEVEL == 2) ~ 
        sum(ifelse(PADLEVEL == 2, PADTIMES * PADDURAT, 0), na.rm = TRUE)/(30/7),
      TRUE ~ NA_real_  # no IAF rows for vigorous
    ),
    mins_moderate_act_wk_earlywaves = case_when(
      all(PAD320 == 0) ~ 0,  # PAQ says no moderate
      any(PADLEVEL == 1) ~ 
        sum(ifelse(PADLEVEL == 1, PADTIMES * PADDURAT, 0), na.rm = TRUE)/(30/7),
      TRUE ~ NA_real_  # no IAF rows for moderate
    ),
    WAVE=WAVE.x,
    .groups = "drop"
  ) %>% group_by(SEQN) %>% slice(1)
# Append to filelist
filelist <- append(filelist, list(new_df))

# Also process dietary/nutrition data separately (see `nhanes-diet-script.R`)
diet_file <- readRDS("nhanes_diet_data/diet_file_final_big.rds")
diet_file$WAVE = factor( diet_file$WAVE )
filelist <- append(filelist, list(diet_file))

# Put all files together
dat <- full_join_reduce_coalesce_clean(filelist, by = "SEQN")
gc()
dat <- dat %>% group_by(SEQN) %>% slice(1) # resolve issue w/ multiples of same row
columnlabels <- column_labels(dat, default = "name")
columnlabels[which(names(columnlabels)=='WAVE')] <- "Wave of Data Collection"



################################################################################
#  Clean outcome variables
################################################################################
# Start with physical activity variables
zeros_PAQ650 <- dat$PAQ650==2
dat$PAQ655[zeros_PAQ650] <- 0;
dat$PAQ655[dat$PAQ655==77] <- NA; 
dat$PAQ655[dat$PAQ655==99] <- NA;
dat$PAD660[zeros_PAQ650] <- 0;
dat$PAD660[dat$PAD660==7777] <- NA;
dat$PAD660[dat$PAD660==9999] <- NA;
zeros_PAQ665 <- dat$PAQ665==2
dat$PAQ670[zeros_PAQ665] <- 0;
dat$PAQ670[dat$PAQ670==77] <- NA;
dat$PAQ670[dat$PAQ670==99] <- NA;
dat$PAD675[zeros_PAQ665] <- 0;
dat$PAD675[dat$PAD675==7777] <- NA;
dat$PAD675[dat$PAD675==9999] <- NA;
dat <- dat %>% mutate(minperwk_vigorousrec=prod(PAQ655,PAD660,na.rm=T),
                      minperwk_moderaterec=prod(PAQ670,PAD675,na.rm=T))
dat$minperwk_vigorousrec[is.na(dat$PAQ655)|is.na(dat$PAD660)] <- NA
dat$minperwk_moderaterec[is.na(dat$PAQ670)|is.na(dat$PAD675)] <- NA
dat <- dat %>%
  mutate(
    minperwk_moderaterec_num  = suppressWarnings(as.numeric(minperwk_moderaterec)),
    minperwk_vigorousrec_num  = suppressWarnings(as.numeric(minperwk_vigorousrec)),
    moderate_leisure = dplyr::coalesce(minperwk_moderaterec_num,
                                       mins_moderate_act_wk_earlywaves),
    vigorous_leisure = dplyr::coalesce(minperwk_vigorousrec_num,
                                       mins_vigorous_act_wk_earlywaves)
  )
# Blood pressure
dat <- dat %>% mutate(sys_bp = mean(c(BPXSY1,BPXSY2,BPXSY3,BPXSY4),na.rm=T))
# Cholesterol
dat <- dat %>% mutate(non_HDL_chol_new = LBDTCSI - LBDHDDSI)
dat <- dat %>% mutate(non_HDL_chol_old = LBDTCSI - LBDHDLSI)
dat <- dat %>% mutate(non_HDL_chol = coalesce(non_HDL_chol_new,
                                              non_HDL_chol_old))
# Body composition
dat <- dat %>% mutate(waist_height = BMXWAIST / BMXHT)
# Creatinine and eGFR
dat <- dat %>% mutate(creatinine = coalesce(LBXSCR,
                                            LBDSCR))
dat <- dat %>% mutate(
  kappa_eGFR = ifelse(RIAGENDR==1,0.9,0.7), #men=1, women=2
  alpha_eGFR = ifelse(RIAGENDR==1,-0.302,-0.241), #men=1, women=2
  female_factor = ifelse(RIAGENDR==1,1,1.012), #men=1, women=2
  # eGFR = 142 x min(Scr/κ, 1)^α x max(Scr/κ, 1)^-1.200 x 0.9938^Age x 1.012[if female]
  eGFR = 142*
    min(creatinine/kappa_eGFR,1)^alpha_eGFR*
    max(creatinine/kappa_eGFR,1)^-1.200*
    0.9938^RIDAGEYR*
    female_factor
) #https://www.frontiersin.org/journals/cardiovascular-medicine/articles/10.3389/fcvm.2024.1417926/full sec 2.2
# changed to https://www.kidney.org/ckd-epi-creatinine-equation-2021 (note removal of race)
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
AGE1_outcome_names <- c('BMXBMI',
                        'waist_height',
                        'LBXCOT',
                        'sys_bp',
                        'non_HDL_chol',
                        diet_outcomes_names)



################################################################################
#  Get and clean confounders, run matching
################################################################################
# treatment is "ratio of family income to poverty" (FIPR)
trtcolumn <- match('INDFMPIR',colnames(dat))
columnlabels[trtcolumn]
dat <- dat[!is.na(dat[trtcolumn]),]
# subset appropriate ages
dat <- dat %>%
  mutate(
    AGEGRP = case_when(
      RIDAGEYR >= 8  & RIDAGEYR <= 11 ~ "1",
      RIDAGEYR >= 12 & RIDAGEYR <= 17 ~ "2",
      TRUE                            ~ NA_character_
    ),
    AGEGRP = factor(AGEGRP, levels = 1:2)
  )
columnlabels[colnames(dat)=='AGEGRP'] <- "Age group of participant"
dat <- dat[dat$RIDAGEYR >= 8 & dat$RIDAGEYR <= 17,]
# get insurance status
dat <- dat %>%
  mutate(
    INSUR = dplyr::coalesce(HID010,HIQ011),
    INSURANCE = case_when(
      INSUR == 1 ~ "1",
      INSUR == 2 ~ "0",
      TRUE           ~ NA_character_
    ),
    INSURANCE = factor(INSURANCE, levels = 0:1)
  )
columnlabels[colnames(dat)=='INSURANCE'] <- "Covered by health insurance?"
dat <- dat %>% filter(!is.na(INSURANCE)) # doesn't actually get rid of many ppl
# only include subjects who were both interviewed and examined
dat <- dat[which(dat$RIDSTATR==2),]
# exclude diagnosed diabetes, currently pregnant (in accordance w previous paper)
dat <- dat %>% filter(DIQ010 == 2)
dat <- dat %>% filter(is.na(RIDEXPRG) | RIDEXPRG==2 )
gc()
# match on gender, race, age, etc.
matchcolumns <- match(c('WAVE','RIDEXMON', 'RIAGENDR',
                        'RIDAGEYR', 'RIDRETH1', 'AGEGRP', 'INSURANCE'),
                      colnames(dat))
columnlabels[matchcolumns]

# define main matching function
run_matching <- function(fipr_cutoff=1, link='probit', cal=0.25, files=filelist){
  # treatment is "Ratio of family income to poverty"
  incometopov <- dat[trtcolumn]
  treatment <- incometopov < fipr_cutoff
  control <- incometopov > fipr_cutoff
  trtvec <- rep(NA,length(incometopov))
  trtvec[control] <- 0
  trtvec[treatment] <- 1
  dat <- dat[!is.na(trtvec),]
  trtvec <- trtvec[!is.na(trtvec)]

  row_missing_outcome <- function(datarow){
    # return boolean, where T means the row is missing any outcome
    if(as.numeric(datarow[names(datarow) %in% c('AGEGRP')])==1){
      datarow=as.numeric(datarow[names(datarow) %in% AGE1_outcome_names])
    }
    any(is.na(datarow))
  }
  
  ixremove <- apply(subset(dat,select=c(all_outcome_names,'AGEGRP')),1,
                    row_missing_outcome) # remove units with any missing outcomes
  dat <- dat[!ixremove,]
  trtvec <- trtvec[!ixremove]
  
  # match on gender, race, age, etc.
  covariates <- dat[,matchcolumns]
  # note factor variables
  covariates[,c(2,3,5)] <- lapply(covariates[,c(2,3,5)],as.factor)
  colnames(covariates) <- c('NHANESWave','TimePd', 'Gender', 'Age',
                            'Ethnicity', 'AgeGrp', 'Insurance')
  df=cbind(covariates,trtvec)
  
  data_matching=covariates
  treated=trtvec
  dff=cbind(data_matching,treated)
  
  datatemp=dff
  datatemp$treated=as.numeric(as.character(dff$treated))
  propscore.model = glm(treated~., family=binomial(link=link), data=datatemp)
  dmy=dummyVars('treated~.',data=datatemp)
  logit.propscore=as.numeric(predict(propscore.model,newx=as.matrix(datatemp)))
  datatemp=cbind(datatemp,logit.propscore)
  
  m.out <- matchit(treated ~ ., exact = ~ NHANESWave+Gender+AgeGrp+Ethnicity,
                   data = datatemp, distance = "robust_mahalanobis",
                   caliper = c(logit.propscore = cal),
                   std.caliper = c(TRUE))
  
  res.m.out <- summary(m.out)
  standardized_difs <- cbind(abs(res.m.out$sum.all[,'Std. Mean Diff.']),
                             abs(res.m.out$sum.matched[,'Std. Mean Diff.']))
  colnames(standardized_difs) <- c('ASD.BF', 'ASD.AF')
  
  treated_idx <- rownames(m.out$match.matrix)
  control_idx <- as.character(m.out$match.matrix)
  treated_rm <- names(m.out$match.matrix[which(is.na(m.out$match.matrix)),])
  treated_idx <- treated_idx[!(treated_idx %in% treated_rm)]
  control_idx <- control_idx[!is.na(control_idx)]
  
  rownames(standardized_difs)[1:(nrow(standardized_difs)-1)] <- 
    c('Wave 99-00','Wave 01-02','Wave 03-04','Wave 05-06','Wave 07-08',
      'Wave 09-10','Wave 11-12','Wave 13-14','Wave 15-16',
      'Exam Nov-Apr','Exam May-Oct',
      'Male','Female',
      'Age',
      'Mex American','Other Hisp','White','Black','Other',
      'Age 8-11','Age 12-17',
      'Uninsured', 'Insured')
  plot.dataframe=data.frame(SMD=c(abs(res.m.out$sum.matched[1:(nrow(res.m.out$sum.matched)-1),
                                                            'Std. Mean Diff.']),
                                  abs(res.m.out$sum.all[1:(nrow(res.m.out$sum.all)-1),
                                                        'Std. Mean Diff.'])),
                            Covariates=rep((rownames(standardized_difs))[1:(nrow(standardized_difs)-1)],2),
                            type=
                              c(rep("After Matching",
                                    length((rownames(standardized_difs))[1:(nrow(standardized_difs)-1)])),
                                rep("Before Matching",
                                    length((rownames(standardized_difs))[1:(nrow(standardized_difs)-1)]))))
  cov_order=(rownames(standardized_difs))[1:(nrow(standardized_difs)-1)]
  plt=ggplot(plot.dataframe,aes(x=SMD,y=Covariates))+
    geom_point(size=3,aes(shape=factor(type,levels = c('Before Matching','After Matching'))))+
    scale_shape_manual(values =c(21,16))+
    scale_y_discrete(limits = rev(cov_order))+
    geom_vline(xintercept=c(0,0.1),lty=2) +
    xlim(0,.4) +
    labs(x = "Absolute Std. Mean Differences", y="Covariates") + theme_bw() +
    theme(axis.text.y=element_text(size=11),
          axis.text.x=element_text(size=11),
          axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11),
          legend.text=element_text(size=11),
          legend.title = element_blank(),legend.position="bottom")
  
  return(list(dat=dat,
              treated_idx=treated_idx,
              control_idx=control_idx,
              num_matched_sets=length(treated_idx),
              bal_tbl=standardized_difs,
              love_plt=plt))
}

# get results
set.seed(0) # for reproducibility
res=run_matching()
res
res$love_plt
dat <- res$dat
treated_idx <- as.numeric(res$treated_idx)
control_idx <- as.numeric(res$control_idx)
num_matched_sets <- res$num_matched_sets


