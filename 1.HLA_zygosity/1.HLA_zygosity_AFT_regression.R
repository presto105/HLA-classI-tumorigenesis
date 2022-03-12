library(dplyr)
library(survival)
library(stringr)

#######################
## defince functions ##
#######################

## define Accelerated failure time model (AFT model)
expected_time_weibull <- function(beta, log_scale){
  return(exp(beta)*gamma(1+exp(log_scale)))}

AFT_regression <- function(data){
  data$tumor <- as.factor(as.character(data$tumor))
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2+tumor"))
  
  m1 <- AIC(survreg(formula=form, data=data, dist="weibull"))
  m2 <- AIC(survreg(formula=form, data=data, dist="exponential"))
  m3 <- AIC(survreg(formula=form, data=data, dist="lognormal"))
  m4 <- AIC(survreg(formula=form, data=data, dist="loglogistic"))
  
  dist_i <- which.min(c(m1, m2, m3, m4))
  print(c(m1, m2, m3, m4))
  print(dist_list[dist_i])
  
  model <- summary(survreg(formula=form, data=data, dist=dist_list[dist_i]))$table
  
  for(i in 1:nrow(model)){
    parameter <- rownames(model)[i]
    time_ratio <- exp(model[i,1])
    lower <- exp(model[i,1]-qnorm(0.975)*model[i,2])
    upper <- exp(model[i,1]+qnorm(0.975)*model[i,2])
    hetero_age <- round(expected_time_weibull(model[1,1], model[length(rownames(model)),1]),2)
    homo_age <- round(expected_time_weibull(model[1,1]+model[i,1], model[length(rownames(model)),1]),2)
    diff_age <- hetero_age-homo_age
    p_value <- model[i,4]
    
    time_ratio <- sprintf("%.3f", time_ratio)
    lower <- sprintf("%.3f", lower)
    upper <- sprintf("%.3f", upper)
    
    time_ratio_with_CI <- paste0(time_ratio," (",lower,",",upper,")")
    p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))
    
    tmp <- data.frame(parameter, time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age)
    
    if(i==1){
      multi_res <- tmp
    }else{
      multi_res <- rbind(multi_res, tmp)
    }
  }
  return(multi_res)
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
         tumor = as.factor(as.character(project_id))
         ) %>% 
  filter(tumor_stage %in% stage1&prior_malignancy=="no"&!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count)) %>% 
  filter(Germline_variant==0&(ILC_Infection_status=="X"|is.na(ILC_Infection_status))) 

## Filter tumors for less than 10 patients
Homo_cnt_over10 <- df %>% group_by(tumor, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n>=10)
tumor_types <- sort(as.character(unique(Homo_cnt_over10$tumor)))

n_by_types <- table(df$tumor)[table(df$tumor)!= 0]
print(n_by_types)

## Calculate AIC values to select individual tumor distributions
for(i in 1:length(tumor_types)){
  tumor <- tumor_types[i]
  AFT_data <- df[which(df$tumor==tumor_types[i]),]
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2"))
  
  m1 <- AIC(survreg(formula=form, data=AFT_data, dist="weibull"))
  m2 <- AIC(survreg(formula=form, data=AFT_data, dist="exponential"))
  m3 <- AIC(survreg(formula=form, data=AFT_data, dist="lognormal"))
  m4 <- AIC(survreg(formula=form, data=AFT_data, dist="loglogistic"))
  
  tmp <- data.frame(tumor, weibull=m1, exponential=m2, lognormal=m3, loglogistic=m4)
  if(i==1){
    res <- tmp 
  }else{
    res <- rbind(res, tmp)
  }
}
print(res)

which_dist <- apply(res[,2:5],1,which.min)
dist_list <- c("weibull","exponential","lognormal","loglogistic")

## AFT regression in the individual tumors with AIC results
for(i in 1:length(tumor_types)){
  AFT_data <- df[which(df$tumor==tumor_types[i]),]
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2"))
  
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist_list[which_dist[i]]))$table

  homo_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='1'))
  hetero_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='0'))
  
  tumor <- tumor_types[i]
  time_ratio <- exp(model[2,1])
  lower <- exp(model[2,1]-qnorm(0.975)*model[2,2])
  upper <- exp(model[2,1]+qnorm(0.975)*model[2,2])
  hetero_age <- round(expected_time_weibull(model[1,1], model[length(rownames(model)),1]),2)
  homo_age <- round(expected_time_weibull(model[1,1]+model[2,1], model[length(rownames(model)),1]),2)
  diff_age <- hetero_age-homo_age
    
  p_value <- model[2,4]
  
  time_ratio <- sprintf("%.3f", time_ratio)
  lower <- sprintf("%.3f", lower)
  upper <- sprintf("%.3f", upper)
  
  time_ratio_with_CI <- paste0(time_ratio," (",lower,",",upper,")")
  p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))

  tmp <- data.frame(tumor, homo_n, hetero_n, time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age)
  if(i==1){
    uni_tumor_res <- tmp
  }else{
    uni_tumor_res <- rbind(uni_tumor_res, tmp)
  }
}
uni_tumor_res <- uni_tumor_res %>% arrange(desc(time_ratio_with_CI))
print(uni_tumor_res)

## AFT regression in the pan-cancer
## The influence of HLA zygosity adjusted with Race and gender
Pan_cancer_df <- df %>% filter(tumor %in% tumor_types)
Pan_cancer_res <- AFT_regression(Pan_cancer_df)
Pan_cancer_res[2,1] <- "Pan_cancer_HLA"

print(Pan_cancer_res)

## AFT regression in the pan-cancer without ccRCC
## The influence of HLA zygosity adjusted with Race and gender
tumor_types <- as.character(tumor_types[!tumor_types %in% 'KIRC'])

Pan_cancer_wo_ccrcc_types_df <- df %>% filter(tumor %in% tumor_types)
Pan_cancer_wo_ccrcc_res <- AFT_regression(data=Pan_cancer_wo_ccrcc_types_df)
Pan_cancer_wo_ccrcc_res[2,1] <- "Pan_cancer_wo_ccRCC_HLA"

print(Pan_cancer_wo_ccrcc_res)

