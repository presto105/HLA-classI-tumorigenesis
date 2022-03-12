library(dplyr)
library(reshape2)
library(ggplot2)
library(survival)

#######################
## defince functions ##
#######################

## define Accelerated failure time model (AFT model) for multiple tumors
expected_time_weibull <- function(beta, log_scale){
  return(exp(beta)*gamma(1+exp(log_scale)))}

dist_list <- c("weibull","exponential","lognormal","loglogistic")
Age_diff_CI_cal <- function(data, form){
  data$tumor <- as.factor(as.character(data$tumor))

  form <- formula(form)
  
  m1 <- AIC(survreg(formula=form, data=data, dist="weibull"))
  m2 <- AIC(survreg(formula=form, data=data, dist="exponential"))
  m3 <- AIC(survreg(formula=form, data=data, dist="lognormal"))
  m4 <- AIC(survreg(formula=form, data=data, dist="loglogistic"))
  
  dist_i <- which.min(c(m1, m2, m3, m4))
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
      age_index_min <<-  floor(min(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
      age_index_max <<- ceiling(max(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
      interval <<- age_index_max-age_index_min
      
      multi_res <- tmp
    }else{
      multi_res <- rbind(multi_res, tmp)
    }
  }
  return(multi_res)
}

## define a visualization function of the count distribution of genotype
genotype_dist <- function(gene, title){
  HLA_genotype <- unique(sapply(strsplit(as.vector(unlist(filtered_df[gene])),','),'[',1))
  HLA_genotype <- sort(HLA_genotype)
  
  genotype_df <- data.frame("genotype" = HLA_genotype)
  for(g in HLA_genotype){
    cnt <- sum(grepl(g, as.vector(unlist(filtered_df[gene]))))
    if (g == HLA_genotype[1]){
      cnt_res <- data.frame("genotype" = g, "count"=cnt)
    }else {
      cnt_res <- rbind(cnt_res, data.frame("genotype" = g, "count"=cnt))
    }
  }
  
  cnt_bar <- cnt_res %>%
    ggplot(aes(x=reorder(genotype,count), y=count)) +
    geom_bar(stat="identity", fill="steelblue")+
    coord_flip() +
    xlab("Allele") + ylab("Allele count") +
    theme_minimal()
  return(cnt_bar)
}

##################
## data process ##
##################

## data load
df_FILE="TCGA_ann_parsed_for_allele.tsv"
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

filtered_df <- df %>% filter(tumor%in%tumor_types)

#########################################
## analysis tumor onset for HLA allele ##
#########################################

## visualization of number distribution of HLA genotypes
genotype_dist('HLA.A', 'HLA-A')
genotype_dist('HLA.B', 'HLA-B')
genotype_dist('HLA.C', 'HLA-C')

## AFT regression in the pan-cancer
ale_status_col <- colnames(filtered_df)[82:157]
alleles=paste0(ale_status_col, collapse='+')
form <- paste0("Surv(age_at_diagnosis,event)~gender+race2+tumor+", alleles, collapse='+')

Pan_cancer_ale_res <- Age_diff_CI_cal(filtered_df, form)
print(Pan_cancer_ale_res)

## AFT regression in the pan-cancer without ccRCC
tumor_types <- as.character(tumor_types[!tumor_types %in% 'KIRC'])
Pan_cancer_wo_ccrcc_types_df <- df %>% filter(tumor %in% tumor_types)

Pan_cancer_wo_ccrcc_ale_res <- Age_diff_CI_cal(data=Pan_cancer_wo_ccrcc_types_df, form)
print(Pan_cancer_wo_ccrcc_ale_res)

## AFT regression in ccRCC
ccrcc_df <- df %>% filter(tumor =='KIRC')
form <- paste0("Surv(age_at_diagnosis,event)~gender+race2+", alleles, collapse='+')

ccrcc_ale_res <- Age_diff_CI_cal(data=ccrcc_df, form)
print(ccrcc_ale_res)

####################################################
## analysis tumor onset for zygosity of HLA genes ##
####################################################

## Visualization of the count distribution of zygosity by HLA gene
counts_df <- melt(rbind(
  `HLA A`=table(filtered_df$A_zygosity),
  `HLA B`=table(filtered_df$B_zygosity),
  `HLA C`=table(filtered_df$C_zygosity)
  ))
colnames(counts_df) <- c('HLA gene', 'Zygosity', 'Counts')

zygosity_dist <- counts_df %>%
  ggplot(aes(x=`HLA gene`, y=Counts, fill=Zygosity)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(text = element_text(size=15))
print(zygosity_dist)

## AFT regression in the pan-cancer
form <- "Surv(age_at_diagnosis,event)~A_zygosity+B_zygosity+C_zygosity+gender+race2+tumor"
Pan_cancer_res <- Age_diff_CI_cal(filtered_df, form)
print(Pan_cancer_res)
write.table(Pan_cancer_res,  "HLA_pan_loglogistic.tsv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

## AFT regression in the pan-cancer without ccRCC
## The influence of individual HLA genes
tumor_types <- as.character(tumor_types[!tumor_types %in% 'KIRC'])

Pan_cancer_wo_ccrcc_types_df <- df %>% filter(tumor %in% tumor_types)
Pan_cancer_wo_ccrcc_res <- Age_diff_CI_cal(data=Pan_cancer_wo_ccrcc_types_df, form)
print(Pan_cancer_wo_ccrcc_res)
write.table(Pan_cancer_wo_ccrcc_res,  "HLA_pan_woccRCC_loglogistic.tsv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


## AFT regression in ccRCC
## The influence of individual HLA genes
ccrcc_df <- df %>% filter(tumor =='KIRC')
form <- "Surv(age_at_diagnosis,event)~A_zygosity+B_zygosity+C_zygosity+gender+race2"
ccrcc_res <- Age_diff_CI_cal(data=ccrcc_df, form)
print(ccrcc_res)

write.table(ccrcc_res,  "HLA_ccRCC_lognormal.tsv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

