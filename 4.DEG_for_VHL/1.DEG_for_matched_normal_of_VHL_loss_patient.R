library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)

#######################
## defince functions ##
#######################

## define a visualization function of GO barplot
Drawing_GO_plot <- function(df, count ,title){
  df_top <- df@result[0:count,]
  GO_barplot <- df_top %>% arrange(desc(p.adjust)) %>% 
    ggplot(aes(x=reorder(Description,-p.adjust),y=Count,fill=p.adjust))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
    geom_bar(stat="identity") +
    xlab("") + ylab("") +
    scale_fill_gradient(low = "#e74c3c", high = "#fadbd8") +
    coord_flip() +
    theme(text = element_text(size=20)) +
    ggtitle(title)
  
  print(GO_barplot)
}

##################
## data process ##
##################

## data load
df_FILE="TCGA_ann.tsv"
raw_df <- read.csv(df_FILE, header=T,sep='\t')

names(raw_df)[names(raw_df)=="age_at_initial_pathologic_diagnosis"] <- "age_at_diagnosis"

## Filter by the criteria described in the method
stage1 <- c("Stage I", "Stage IB", "Stage IA", "Stage IC", "Stage IB1", "Stage 0", "Stage IS", "Stage IB2")
df <- raw_df %>%
  mutate(event = 1,
         Homo_at_least_1_locus = as.factor(ifelse(Homozygous_gene_count>=1,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Germline_variant = as.factor(ifelse(P_total_GV_count==0,0,1)),
         age_at_diagnosis = as.numeric(as.character(age_at_diagnosis)),
         prior_malignancy = as.factor(as.character(prior_malignancy)),
         tumor = as.factor(as.character(project_id)),
         VHL = relevel(as.factor(ifelse(is.na(VHL_mut_type), NA,
                                        ifelse(VHL_mut_type=='None','None','Mutation'))),ref='None'),
         chr3_VHL = relevel(as.factor(ifelse(is.na(chr3p_VHL_binned), NA,
                                             ifelse(chr3p_VHL_binned=='-1','-1','0'))),ref='0'),
         VHL_biallelic = factor(ifelse(VHL=='Mutation'&chr3p_VHL_binned==-1,'loss','None'),levels=c('None','loss')),
         VHL_loss = factor(ifelse(VHL=='Mutation'|chr3p_VHL_binned==-1,'loss','None'),levels=c('None','loss'))
  ) %>% 
  filter(tumor_stage %in% stage1&prior_malignancy=="no"&!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count)) %>% 
  filter(Germline_variant==0&(ILC_Infection_status=="X"|is.na(ILC_Infection_status))) %>% 
  filter(tumor=='KIRC'&!is.na(VHL_biallelic))

########## Data preprocessing ##########
KIRC_exp_FILE="KIRC.mRNAseq_raw_counts_tumors.txt"
KIRC_exp_df <- read.csv(KIRC_exp_FILE, header = TRUE, sep='\t')

## selecting samples that exist a matched normal
VHL_loss_pts <- (df %>% filter(VHL_biallelic=='loss'))$submitter_id
VHL_loss_pts <- gsub("-",".",VHL_loss_pts)
VHL_loss_pts <- colnames(KIRC_exp_df)[substring(colnames(KIRC_exp_df), 1, 12)%in% VHL_loss_pts]

VHL_loss_targets_df <- KIRC_exp_df[c('gene_id', VHL_loss_pts)]
matched_normal_ids <- colnames(VHL_loss_targets_df)[substr(colnames(VHL_loss_targets_df),14,15)=="11"] ## select normal sampels
matched_samples <- colnames(VHL_loss_targets_df)[substr(colnames(VHL_loss_targets_df),1, 12) %in% substr(matched_normal_ids, 1, 12)]

## parsing dataframe for gene count
VHL_NT_df <- VHL_loss_targets_df[c('gene_id', matched_samples)]
rownames(VHL_NT_df) <- VHL_NT_df$gene_id
VHL_NT_df <- VHL_NT_df[,-1]

## generating data annotated as normal and tumor
coldata <- data.frame("sample" = colnames(VHL_NT_df))
coldata <- coldata %>%
  mutate(
    condition = ifelse(substr(sample, 14,15)=='11','normal','tumor')
    )
VHL_mat <- as.matrix(VHL_NT_df)

##################
## RNA analysis ##
##################

## running DeSeq2 for DEG anlaysis
dds <- DESeqDataSetFromMatrix(countData = VHL_mat,
                            colData = coldata,
                            design=~condition)
dds <- DESeq(dds)
res <- results(dds)
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
keep <- rowSums(nc>=1) >= (ncol(VHL_mat)/2)
dds<- dds[keep,]
dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
mat <- assay(vst)
res <- results(dds)

res$padj[is.na(res$padj)] <- 1
filtered_res <- res[res$padj <= 0.05, ]

print(filtered_res)

## pathway analysis
LFC_list <- list()
geneList <- list()
up_down_mat_names <- c("up", "down")

target_genes <- rownames(res)
target_genes <- sapply(strsplit(target_genes, "\\|"), "[[", 1)
target_genes <- target_genes[!duplicated(target_genes)]

geneList_LFC <- res$log2FoldChange
names(geneList_LFC) <- as.character(target_genes)
geneList_LFC <- sort(geneList_LFC, decreasing = TRUE)
LFC_list[[1]] <- geneList_LFC[geneList_LFC > 1]
LFC_list[[2]] <- geneList_LFC[geneList_LFC < -1]
geneList[[1]] <- names(geneList_LFC)[geneList_LFC > 1]
geneList[[2]] <- names(geneList_LFC)[geneList_LFC < -1]

## translate biological ID
ENS.geneList <- list()
for (k in 1:length(geneList)){
  ENS.geneList[[k]] <- bitr(geneList[[k]], fromType = "SYMBOL",
                                          toType = c("ENSEMBL", "ENTREZID"),
                                          OrgDb = org.Hs.eg.db)
  ENS.geneList[[k]] <- ENS.geneList[[k]][,c(3,2,1)]
  
}  

## GO enrichment analysis of a gene set
ego_comp <- list()
ego_filter_comp <- list()
for (y in 1:length(ENS.geneList)) {
    ego_comp[[y]] <- enrichGO(gene = ENS.geneList[[y]]$SYMBOL,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType = "SYMBOL",
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05
  )
}

Drawing_GO_plot(ego_comp[[1]], 15, "GENE ONTOLOGY UP")
Drawing_GO_plot(ego_comp[[2]], 15, "GENE ONTOLOGY DOWN")
