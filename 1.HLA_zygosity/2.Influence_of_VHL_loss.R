library(dplyr)
library(survival)

#######################
## defince functions ##
#######################

## define Accelerated failure time model (AFT model)
expected_time_weibull <- function(beta, log_scale){
  return(exp(beta)*gamma(1+exp(log_scale)))}

dist_select <- function(data, var){
  AFT_data <- data
  form <- formula(paste("Surv(age_at_diagnosis,event)~",paste(var,collapse="+")))

  AIC_res <- c()
  dist_class <- c("weibull", "exponential", "lognormal", "loglogistic")
  for(dist in dist_class){
    AIC_res <- c(AIC_res, AIC(survreg(formula=form, data=AFT_data, dist=dist)))
  }
  AIC_df <- data.frame(AIC_res, row.names=dist_class)
  min_dist <- dist_class[which(AIC_res==min(AIC_res))]
  
  return(list(min_dist=min_dist,AIC_df=AIC_df))
}

AFT_regression <- function(data, dist){
  AFT_data <- data

  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2+tumor"))
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))$table
  
  n <- nrow(AFT_data)
  homo_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='1'))
  hetero_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='0'))
  
  Time_ratio <- exp(model[2,1])
  lower <- exp(model[2,1]-qnorm(0.975)*model[2,2])
  upper <- exp(model[2,1]+qnorm(0.975)*model[2,2])

  hetero_age <- round(expected_time_weibull(model[1,1], model[3,1]),2)
  homo_age <- round(expected_time_weibull(model[1,1]+model[2,1], model[3,1]),2)
  diff_age <- hetero_age-homo_age
  
  p_value <- model[2,4]
  
  Time_ratio <- sprintf("%.3f", Time_ratio)
  lower <- sprintf("%.3f", lower)
  upper <- sprintf("%.3f", upper)
  
  Time_ratio_with_CI <- paste0(Time_ratio," (",lower,",",upper,")")
  p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))
  
  res <- data.frame(homo_n, hetero_n, Time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age)
  return(res)
}

##################
## data process ##
##################

## data load
df_FILE="TCGA_ann.tsv"
raw_df <- read.csv(df_FILE, header=T,sep='\t')

names(raw_df)[names(raw_df)=="age_at_initial_pathologic_diagnosis"] <- "age_at_diagnosis"

## filter by the criteria described in the method
stage1 <- c("Stage I", "Stage IB", "Stage IA", "Stage IC", "Stage IB1", "Stage 0", "Stage IS", "Stage IB2")

df <- raw_df %>% 
  mutate(event = 1,
         Homo_at_least_1_locus = as.factor(ifelse(Homozygous_gene_count>=1,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Germline_variant = as.factor(ifelse(P_total_GV_count==0,0,1)),
         
         age_at_diagnosis = as.numeric(as.character(age_at_diagnosis)),
         race2 = as.factor(ifelse(race=="WHITE",1,0)),
         gender = as.factor(as.character(gender)),
         prior_malignancy = as.factor(as.character(prior_malignancy)),
         tumor = as.factor(as.character(project_id)),
         VHL = relevel(as.factor(ifelse(is.na(VHL_mut_type), NA,
                                        ifelse(VHL_mut_type=='None','None','Mutation'))),ref='None'),
         PBRM1 = relevel(as.factor(ifelse(is.na(PBRM1_mut_type),NA,
                                          ifelse(PBRM1_mut_type=='None','None','Mutation'))),ref='None'),
         chr3p_VHL_binned = relevel(as.factor(ifelse(is.na(chr3p_VHL_binned), NA,
                                                     ifelse(chr3p_VHL_binned=='-1','-1','0'))),ref='0'),
         chr3p_PBRM1_binned = relevel(as.factor(ifelse(is.na(chr3p_PBRM1_binned), NA,
                                                       ifelse(chr3p_PBRM1_binned=='-1','-1','0'))),ref='0'),
         VHL_biallelic = factor(ifelse(VHL=='Mutation'&chr3p_PBRM1_binned==-1,'O','X'),levels=c('X','O'))
  ) %>% 
  filter(prior_malignancy=="no"&!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count)) %>% 
  filter(Germline_variant==0&(ILC_Infection_status=="X"|is.na(ILC_Infection_status))) %>% 
  filter(VHL_vaf>0.1|is.na(VHL_vaf))

###########################################################
## Analysis of tumor onset time in the total tumor stage ##
###########################################################
Homo_cnt_over10 <- df %>% group_by(tumor, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n>=10)
tumor_types <- sort(as.character(unique(Homo_cnt_over10$tumor)))

TS_df <- df %>%
  filter(tumor %in% tumor_types)

## AFT regression of ccRCC (total stage)
TS_ccRCC_df <- TS_df %>% filter(tumor=='KIRC')
dist <- dist_select(TS_ccRCC_df, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_ccRCC_df, dist=dist)

## AFT regression of pan-cancer without ccRCC (total stage)
TS_wo_ccRCC <- TS_df %>% filter(tumor!='KIRC')
dist <- dist_select(TS_wo_ccRCC, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_wo_ccRCC, dist=dist)

## AFT regression of VHL bi-allelic loss tumors (total stage without limiting tumor type)
TS_VHL_loss <- TS_df %>% filter(VHL_biallelic=='O')
dist <- dist_select(TS_VHL_loss, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_VHL_loss, dist=dist)

## AFT of non VHL bi-allelic tumors (total stage, without limiting tumor type)
TS_non_VHL_loss <- TS_df %>% filter(VHL_biallelic!='O'|is.na(VHL_biallelic))
dist <- dist_select(TS_non_VHL_loss, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_non_VHL_loss, dist=dist)


## AFT of non VHL bi-allelic tumors (total stage, ccRCC)
TS_non_VHL_loss_ccRCC <- TS_df %>%
  filter(tumor=='KIRC') %>% 
  filter(VHL_biallelic!='O'|is.na(VHL_biallelic))
dist <- dist_select(TS_non_VHL_loss_ccRCC, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_non_VHL_loss_ccRCC, dist=dist)


## AFT of VHL bi-allelic tumors (total stage, ccRCC)
TS_VHL_loss_ccRCC <- TS_df %>%
  filter(tumor!='KIRC') %>% 
  filter(VHL_biallelic=='O')
dist <- dist_select(TS_VHL_loss_ccRCC, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=TS_VHL_loss_ccRCC, dist=dist)


#################################################
## Analysis of tumor onset time in the stage 1 ##
#################################################
S1_df <- df %>%  filter(tumor_stage %in% stage1)

Homo_cnt_over10 <- S1_df %>% group_by(tumor, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n>=10)
tumor_types <- sort(as.character(unique(Homo_cnt_over10$tumor)))

S1_df <- S1_df %>% 
  filter(tumor %in% tumor_types)

## AFT regression of ccRCC (stage 1)
S1_ccRCC_df <- S1_df %>% filter(tumor=='KIRC')
dist <- dist_select(S1_ccRCC_df, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_ccRCC_df, dist=dist)

## AFT regression of pan-cancer without ccRCC (stage 1)
S1_VHL_loss<- S1_df %>% filter(tumor!='KIRC')
dist <- dist_select(S1_VHL_loss, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_VHL_loss, dist=dist)


## AFT regression of VHL bi-allelic loss tumors (stage 1 without limiting tumor type)
S1_VHL_loss<- S1_df %>% filter(VHL_biallelic=='O')
dist <- dist_select(S1_VHL_loss, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_VHL_loss, dist=dist)

## AFT of non VHL bi-allelic tumors (stage 1, without limiting tumor type)
S1_non_VHL_loss<- S1_df %>% filter(VHL_biallelic!='O'|is.na(VHL_biallelic))
dist <- dist_select(S1_non_VHL_loss, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_non_VHL_loss, dist=dist)

## AFT of non VHL bi-allelic tumors (stage 1, ccRCC)
S1_non_VHL_loss_ccRCC<- S1_df %>%
  filter(tumor=='KIRC') %>% 
  filter(VHL_biallelic!='O'|is.na(VHL_biallelic))
dist <- dist_select(S1_non_VHL_loss_ccRCC, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_non_VHL_loss_ccRCC, dist=dist)

## AFT of VHL bi-allelic tumors (stage 1, ccRCC)
S1_VHL_loss_ccRCC<- S1_df %>%
  filter(tumor!='KIRC') %>% 
  filter(VHL_biallelic=='O')
dist <- dist_select(S1_VHL_loss_ccRCC, "Homo_at_least_1_locus")[[1]]
AFT_regression(data=S1_VHL_loss_ccRCC, dist=dist)

