library(Rsubread)
library(GEOquery)
library(refGenome)
library(GenomicRanges)
library(Rsamtools)
library(doParallel)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(org.Hs.eg.db)
library(GEOquery)
library(sva)
library(affy)
library(limma)
library(DESeq2)
library(glmnet)
library(TissueEnrich)
library(clusterProfiler)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(survivalROC)
library(randomForest)
library(randomForestSRC)
library(factoextra)
library(ggfortify)
library(ranger)
library(ggRandomForests)
library(maftools)
library(RColorBrewer)
library(gghalves)
library(ggsci)
library(ggalluvial)
library(ggpubr)
library(ggridges)
library(magrittr)
library(furrr)



library(tidyverse)


rm(list = ls())


# 1. RAW_microarray_Limma -------------------------------------------------


raw_dir = c("./input/raw_data/Con1.txt", "./input/raw_data/Con2.txt","./input/raw_data/Con3.txt", 
  "./input/raw_data/Exp1.txt","./input/raw_data/Exp2.txt", "./input/raw_data/Exp3.txt")
raw_data <- read.maimages(raw_dir,source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")

prob_filter <- openxlsx::read.xlsx("./input/raw_data/Matrix.xlsx")

prob_gene <- read.csv(file = "./input/raw_data/probe2ensembl_GPL26963_v36.csv")[,2:3] %>% 
  rename("ProbeName" = "Probe", "GeneID"="Ensembl")

expr_data <- backgroundCorrect(raw_data, method="normexp") %>% 
  normalizeBetweenArrays(method="quantile")
Control <- expr_data$genes$ControlType==1L
IsExpr <- rowSums(expr_data$other$gIsWellAboveBG > 0) >= 3
expr_data <- expr_data[!Control & IsExpr, ]

sample_infor <- data.frame(sample_name = c('Con1','Con2','Con3','Exp1', 'Exp2','Exp3'),
  Batch = as.factor(c('A','B','C','A','B','C')),
  disease = c('Con','Con','Con','Exp','Exp','Exp'))
mod = model.matrix(~as.factor(sample_infor$disease))

expr <- expr_data$E
rownames(expr) <- expr_data[["genes"]][["ProbeName"]]
colnames(expr) <- c('Con1','Con2','Con3', 'Exp1', 'Exp2','Exp3')
expr %<>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ProbeName") %>% 
  filter(ProbeName %in% prob_filter$ID_REF) %>% 
  column_to_rownames(var = "ProbeName") %>% 
  as.data.frame() %>% 
  ComBat(batch = sample_infor$Batch,mod = mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ProbeName") %>% 
  filter(ProbeName %in% prob_gene$ProbeName) %>% 
  left_join(prob_gene[,1:2],by = "ProbeName") %>% 
  column_to_rownames(var = "ProbeName") %>% 
  select("GeneID",everything()) %>% 
  as.data.frame()

expr <- aggregate(x = expr[,2:ncol(expr)],
  by = list(expr$GeneID),
  FUN = max) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "Group.1")

pair <- factor(rep(1:3, each=2))
group <- factor(c(rep("Con",3),rep('Exp',3)))
pair_design <- model.matrix(~pair+group)
colnames(pair_design)

row_limma <- topTable(eBayes(lmFit(expr, pair_design)), coef="groupExp",n = Inf, adjust.method = "fdr") %>% 
  rownames_to_column(var = "GeneID") %>% 
  left_join(prob_gene[,1:2],by = "GeneID") %>% 
  filter(!duplicated(GeneID)) %>% 
  filter((abs(logFC) > log2(1.5)) & (P.Value < 0.05)) %>% 
  select("logFC","P.Value","GeneID") %>% 
  #column_to_rownames(var = "GeneSymbol") %>% 
  arrange(desc(logFC))
gene_high <- row_limma %>% 
  filter(logFC > 0)
gene <- row_limma$GeneID[str_detect(row_limma$GeneID,pattern = 'ENSG')]
gene_id <- gene_high$GeneID[str_detect(gene_high$GeneID,pattern = 'ENSG')]

# 2. Figure2B -------------------------------------------------------------


heatmap_expr <- expr %>% 
  rownames_to_column(var = "GeneID") %>% 
  filter(GeneID %in% all_of(row_limma$GeneID)) %>% 
  left_join(prob_gene[,1:2],by = "GeneID") %>% 
  filter(!duplicated(GeneID)) %>% 
  column_to_rownames(var = "GeneID") %>% 
  select("Con1","Con2","Con3","Exp1","Exp2","Exp3")

MicroArray_heatmap <- pheatmap::pheatmap(heatmap_expr, 
  main = "Heatmap of exosome arry",show_colnames = TRUE,show_rownames=T,
  fontsize=10,color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50),
  annotation_legend=TRUE,border_color=NA,scale="row",cluster_rows = TRUE,cluster_cols = F)
labels_row <- MicroArray_heatmap[["tree_row"]][["labels"]]
labels_row[which(!(labels_row %in% c("ENSG00000065911.13","ENSG00000087842.11",
      "ENSG00000140992.19","ENSG00000176871.9" )))] <- ""

#Figure2B
plot2B <- pheatmap::pheatmap(heatmap_expr, 
  labels_row = labels_row,
  main = "Heatmap of Different Gene Expression in exosomes",
  show_colnames = TRUE,
  fontsize=10,
  color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50),
  annotation_legend=TRUE,
  border_color=NA,
  scale="row",
  treeheight_row = 45,
  cluster_rows = TRUE,
  cluster_cols = F)
plot2B

# 3. Figure2C -------------------------------------------------------------


scatter_exp <- topTable(eBayes(lmFit(expr, pair_design)), coef="groupExp",n = Inf, adjust.method = "fdr") %>% 
  rownames_to_column(var = "GeneID") %>% 
  left_join(prob_gene[,1:2],by = "GeneID") %>% 
  filter(!duplicated(GeneID)) %>% 
  mutate(LogP = -log10(P.Value),
    Group = ifelse(((P.Value > 0.05) | (abs(logFC) < log2(1.5))),"None",
      ifelse((logFC < log2(1.5)) & (P.Value < 0.05),"Low","High")),
    label = GeneID)
library(ggsci)

#Figure2C
plot2C <- ggscatter(scatter_exp,x="logFC", y="LogP",
  color = 'Group',palette = c('#ff0000','#0000ff','#808080'),alpha = 0.6,
  size = abs(scatter_exp$logFC),font.label = 8,repel = T,
  xlim = c(-3,3),ylim = c(0,4))+
  theme_classic2()+
  geom_hline(yintercept=1.30,linetype ="dashed")+
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),linetype ="dashed")+
  labs(title="Different Expression by Array", x="FoldChange", y="-Log10(P-value)")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5)) 
plot2C

# 4. ConsensusClusterPlus & Figure2DEF ------------------------------------


rm(list = setdiff(ls(),c('gene','gene_id')))

#load TCGA data
load("./input/rdata/fpkm_survival.rdata")

fpkm_ConsensusCluster <- fpkm_survival %>% 
  as.data.frame() %>% 
  select(any_of(gene_id)) %>% 
  t()

#Figure2DE
counts_results <- ConsensusClusterPlus(fpkm_ConsensusCluster, maxK = 10, 
  reps = 1000, pItem = 0.8, pFeature = 1, 
  title = "Exosomes_ConsensusCluster", verbose = F, 
  clusterAlg = "km", distance="euclidean", 
  plot = "pdf", writeTable = F)

#proportion of am-biguous clustering (PAC) score shows an optimal K is two
if(T){
  Kvec = 2:10
  x1 = 0.1; x2 = 0.9
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="")
  for(i in Kvec){
    M = counts_results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  optK = Kvec[which.min(PAC)]
  optK
}

consensusClass <- ifelse((counts_results[[2]][["consensusClass"]]) == 1,'exoB','exoA')

exo_class <- cbind(sample_id = names(consensusClass),
  Group = consensusClass) %>% 
  as.data.frame() %>% 
  left_join(data.frame(sample_id = rownames(fpkm_survival),fpkm_survival[,c(1,2)])) %>% 
  column_to_rownames(var = "sample_id") 

TCGA_fit_OS <- survfit(Surv(OS.Time, OS) ~Group, data = exo_class)

#Figure2F
plot2F <- ggsurvplot(TCGA_fit_OS, data = exo_class, 
  risk.table = T, conf.int = F ,pval = T,
  xlab = "Follow up time(DAY)",
  ylab = "OS Probability",
  break.x.by = 365,
  xlim = c(0,3650))
plot2F

# 5. Figure3AB ------------------------------------------------------------


KM_data <- fpkm_survival %>% 
  as.data.frame() %>% 
  select("OS","OS.Time",any_of(gene_id)) 

coxR <- KM_data[,3:ncol(KM_data)] %>% 
  map(function(a){
    mycox <- coxph(Surv(KM_data[,2],KM_data[,1]) ~ as.numeric(a))
    coxResult <- summary(mycox)
    coxR <- cbind(HR = coxResult$coefficients[,"exp(coef)"],
      P = coxResult$coefficients[,"Pr(>|z|)"])
    return(coxR)
  }) %>% 
  reduce(rbind) %>% 
  cbind(colnames(KM_data)[3:ncol(KM_data)]) %>% 
  as.data.frame() %>% 
  dplyr::rename('gene_id' = 'V3')

fpkm_TCGA_gene <- coxR %>% filter(P < 0.05 & as.numeric(HR) > 1)
U133_prob_gene <- read.csv(file = "./input/raw_data/probe2ensembl_HG-U133A_v36.csv")
#we excluded genes that were not available in the Affymetrix Human Genome U133A Array
fpkm_intersect_gene <- intersect(fpkm_TCGA_gene$gene_id,U133_prob_gene$Ensembl)
rf_data <- fpkm_survival %>% 
  as.data.frame() %>% 
  select("OS","OS.Time",all_of(fpkm_intersect_gene))
fpkm_rf_fit <- rfsrc(Surv(OS.Time, OS) ~ ., data = rf_data, ntree = 5000,importance = TRUE)

#Figure3AB
plot(fpkm_rf_fit)

#our genes to construct the final model
fpkm_rf_gene_name <- cbind(gene_id = names(fpkm_rf_fit$importance),score = fpkm_rf_fit[["importance"]]) %>% 
  as.data.frame() %>% 
  filter(score > 0) %>% 
  select(gene_id)


# 6.Figure3CDE ------------------------------------------------------------


rm(list = setdiff(ls(),c('fpkm_rf_gene_name','fpkm_survival')))
fpkm_survival <- fpkm_survival %>% 
  as.data.frame() %>% 
  select("OS","OS.Time",all_of(fpkm_rf_gene_name[,1]))

if(T){
  res_OScox <- coxph(Surv(OS.Time, OS) ~ ., data = fpkm_survival)
  summary(res_OScox)
  
  res <- summary(res_OScox)[["coefficients"]] %>% 
    as.data.frame() 
  
  fpkm_rf_gene <- cbind(str_remove_all(rownames(res),pattern = '[`]'),res$coef) %>% 
    as.data.frame() %>% 
    rename(genename = V1, coef = V2)
}

# ES

tcga_score <- fpkm_survival[,-c(1,2)] %>% 
  select(fpkm_rf_gene[,1]) %>% 
  apply(.,1,function(x){
    tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(fpkm_rf_gene[,2]))
    return(tt)
  })

tcga_clust <- cbind(sample = names(tcga_score), score = as.numeric(tcga_score)) %>% 
  as_tibble() %>% 
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
    OS.Time = as.numeric(fpkm_survival$OS.Time),
    OS = as.numeric(fpkm_survival$OS))

TCGA_fit_OS <- survfit(Surv(OS.Time, OS) ~Group, data = tcga_clust)

#Figure3C
plot3C <- ggsurvplot(TCGA_fit_OS, data = tcga_clust, 
  risk.table = T, conf.int = F,pval = T,
  xlab = "Follow up time(DAY)", ylab = "OS Probability",
  break.x.by = 365, xlim = c(0,4000))
plot3C

multi_ROC <- function(time_vector,time_name,risk_score_table){
  single_ROC <- function(single_time,single_time_name){
    for_ROC <- survivalROC(Stime = risk_score_table$OS.Time,
      status = risk_score_table$OS,
      marker = risk_score_table$score,
      predict.time = single_time, method = 'KM')
    data.frame('Time_point'=rep(single_time_name, length(for_ROC$TP)),
      'True_positive'=for_ROC$TP,
      'False_positive'=for_ROC$FP,
      'Cut_values'=for_ROC$cut.values,
      'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- pmap(list(time_vector, time_name), single_ROC) %>% 
    reduce(rbind) %>% 
    mutate(Time_point = factor(.$Time_point,levels = time_name))
}
time_vector <- 365*c(3,5,7)
time_name <- c('3y','5y','7Y')
TCGA_OS_clust_ROC <- multi_ROC(time_vector = time_vector,time_name = time_name,risk_score_table = tcga_clust)

#Figure3D
plot3D <- ggplot(TCGA_OS_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
    round(TCGA_OS_clust_ROC$AUC[which(TCGA_OS_clust_ROC$Time_point == time_name[1])][1],3)),
    paste("AUC of 5 year =",
      round(TCGA_OS_clust_ROC$AUC[which(TCGA_OS_clust_ROC$Time_point == time_name[2])][1],3)),
    paste("AUC of 7 year =",
      round(TCGA_OS_clust_ROC$AUC[which(TCGA_OS_clust_ROC$Time_point == time_name[3])][1],3)))) +
  theme(strip.background = element_rect(fill="grey90"),
    strip.text = element_text(size=13,face="plain",color="black"),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=11,face="plain",color="black"),
    panel.background=element_rect(colour="black",fill=NA),
    panel.grid=element_blank(), legend.text = element_text(size = 16),
    legend.position=c(0.95,0.05), legend.justification = c(1,0),
    legend.background=element_blank(), legend.key = element_blank(),
    axis.ticks=element_line(colour="black"))
plot3D


load("./input/rdata/tcga_pscore.rdata")

tcga_multiCOX <- tcga_pscore[,1:2] %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(tcga_clust) %>% 
  select("age","stage","score","Group","OS.Time","OS")

TCGA_multiCOX_OS <- tcga_multiCOX %>%
  filter(!(stage == "NA")) %>% 
  filter(!(is.na(OS.Time) | is.na(OS))) %>% 
  filter(!(OS.Time == 0 | OS.Time == '#N/A')) %>% 
  mutate(OS.Time = as.numeric(OS.Time),
    score = as.numeric(score),
    age = as.numeric(age),
    stage = factor(stage,labels = c('I','II','III','IV')),
    OS = as.numeric(OS)) %>% 
  select(OS.Time,OS,age,stage,score)

TCGA_OS_model <- coxph( Surv(OS.Time, OS) ~  age + stage + score , data =  TCGA_multiCOX_OS)

#Figure3E
plot3E <- ggforest(TCGA_OS_model,data = TCGA_multiCOX_OS,
  main = 'Hazard ratio of BRCA (OS)',  cpositions = c(0.05, 0.15, 0.35),  
  fontsize = 1, noDigits = 3)+
  theme_void()+
  theme(strip.background = element_rect(fill="grey90"),
    strip.text = element_text(size=13,face="plain",color="black"),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=11,face="plain",color="black"),
    panel.background=element_rect(colour="black",fill=NA),
    panel.grid=element_blank(),
    legend.position="none",
    legend.background=element_rect(colour=NA,fill=NA),
    axis.ticks=element_line(colour="black"))
plot3E

# 7. Figure3FG ------------------------------------------------------------


load("./input/rdata/GSE25066.rdata")

geo_score <- fpkm_GEO_survival[,-c(1,2)] %>% 
  select(fpkm_rf_gene[,1]) %>% 
  apply(.,1,function(x){
    tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(fpkm_rf_gene[,2]))
    return(tt)
  })

geo_clust <- cbind(sample = names(geo_score), score = as.numeric(geo_score)) %>% 
  as_tibble() %>% 
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
    OS.Time = as.numeric(fpkm_GEO_survival$OS.Time),
    OS = as.numeric(fpkm_GEO_survival$OS))
GEO_fit_OS <- survfit(Surv(OS.Time, OS) ~Group, data = geo_clust)

#Figure3F
plot3F <- ggsurvplot(GEO_fit_OS, data = geo_clust, 
  risk.table = T, conf.int = F,pval = T,
  xlab = "Follow up time(DAY)", ylab = "DRFS Probability",
  break.x.by = 365, xlim = c(0,2920))
plot3F

time_vector <- 365*c(3,5)
time_name <- c('3y','5y')
GEO_OS_clust_ROC <- multi_ROC(time_vector = time_vector,
  time_name = time_name,risk_score_table = geo_clust)

#Figure3G
plot3G <- ggplot(GEO_OS_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
    round(GEO_OS_clust_ROC$AUC
      [which(GEO_OS_clust_ROC$Time_point == time_name[1])][1],3)),
    paste("AUC of 5 year =",
      round(GEO_OS_clust_ROC$AUC
        [which(GEO_OS_clust_ROC$Time_point == time_name[2])][1],3))
  )) +
  theme(strip.background = element_rect(fill="grey90"),
    strip.text = element_text(size=13,face="plain",color="black"),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=11,face="plain",color="black"),
    panel.background=element_rect(colour="black",fill=NA),
    panel.grid=element_blank(),
    legend.text = element_text(size = 16),
    legend.position=c(0.95,0.05),
    legend.justification = c(1,0),
    legend.background=element_blank(),
    legend.key = element_blank(),
    axis.ticks=element_line(colour="black"))
plot3G





# 8. Figure4A -------------------------------------------------------------


rm(list = setdiff(ls(),c('fpkm_rf_gene','fpkm_survival','tcga_pscore')))
Immune_Subtype <- tcga_pscore[,c(3,5)] %>% 
  filter(!(Immune.Subtype == "NA")) %>% 
  mutate(Immune.Subtype = factor(Immune.Subtype),
    ES = as.numeric(score)) %>% 
  select(Immune.Subtype,ES)

my_comparisons <- list( c("C1", "C2"), c("C1", "C3"), c("C1", "C4"), c("C1", "C6"),
  c("C2", "C3"), c("C2", "C4"), c("C2", "C6"),
  c("C3", "C4"), c("C3", "C6"),
  c("C4", "C6"))

#Figure4A
plot4A <- ggplot(Immune_Subtype , aes(x = Immune.Subtype, 
  y = ES, fill = Immune.Subtype))+
  geom_boxplot(aes(x = Immune.Subtype,y = ES, fill = Immune.Subtype),
    outlier.shape = NA,
    width = .07,
    color = "black")+
  stat_compare_means(label.y = 12,size = 5)+
  stat_compare_means(comparisons = my_comparisons)+
  geom_half_violin(aes(fill = Immune.Subtype),
    position = position_nudge(x = .15, y = 0),
    width = .2,
    adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
  geom_point(aes(x = as.numeric(Immune.Subtype)-0.1,
    y = ES,color = Immune.Subtype),
    position = position_jitter(width = .05),size = .25, shape = 20) +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()+
  theme(panel.background=element_rect(colour="black",fill=NA),
    panel.grid=element_blank(),
    legend.position="none",
    legend.background=element_rect(colour=NA,fill=NA))


# 9. Figure4BC -----------------------------------------------------------


immune_factor <- tcga_pscore %>% 
  filter(!is.na(Group)) %>% 
  select(-c("age","stage","Immune.Subtype","Wound.Healing","Macrophage.Regulation",
    "IFN-gamma.Response","TGF-beta.Response",
    "Th1.Cells","Th2.Cells","Th17.Cells","Plasma.Cells","T.Cells.Follicular.Helper",
    "T.Cells.gamma.delta","T.Cells.Regulatory.Tregs",
    "Lymphocyte.Infiltration.Signature.Score","Silent.Mutation.Rate","Nonsilent.Mutation.Rate"))

#data prepare
if(T){
  #cor
  cor_result <- map(immune_factor[,-c(1,2)],function(x){
    #x = immune_factor[,-c(1,2)][,18]
    x = as.numeric(x)
    x = (x - min(x,na.rm = T))/(max(x,na.rm = T) - min(x,na.rm = T))
    y = as.numeric(immune_factor$score)
    y = (y - min(y,na.rm = T))/(max(y,na.rm = T) - min(y,na.rm = T))
    cor_res = cor.test(x,y,exact = F)
    a = cor_res[["p.value"]]
    b = cor_res[["estimate"]][["cor"]]
    re = c(a,b)
    names(re) = c("pvalue","cor")
    return(re)
  }) %>% 
    purrr::reduce(rbind) %>% 
    as.data.frame() 
  
  rownames(cor_result) <- colnames(immune_factor[,c(3:ncol(immune_factor))])
  cor_result$qvalue <- p.adjust(cor_result$pvalue)
  
  
  #T-test
  ttest_result <- map(immune_factor[,-c(1,2)],function(x){
    #x = immune_factor[,-c(1,2)][,18]
    x = as.numeric(x)
    x = (x - min(x,na.rm = T))/(max(x,na.rm = T) - min(x,na.rm = T))
    wilcox_res = t.test(x~immune_factor$Group, exact = F)
    a = wilcox_res[["p.value"]]
    b = wilcox_res[["statistic"]]
    c = cbind(x,immune_factor$Group) %>% 
      as.data.frame() %>% 
      group_by(V2) %>%
      summarize(mean = mean(as.numeric(x),na.rm = T))
    re = c(a,b,c$mean)
    names(re) = c("pvalue","t","mean_High","mean_Low")
    return(re)
  }) %>% 
    purrr::reduce(rbind) %>% 
    as.data.frame() %>% 
    mutate(res = ifelse((as.numeric(.$mean_High) - as.numeric(.$mean_Low))>0,'High','Low')) %>% 
    select("mean_High","mean_Low","t","pvalue","res")
  
  rownames(ttest_result) <- colnames(immune_factor[,c(3:ncol(immune_factor))])
  
  ttest_result %<>% 
    mutate(qvalue = p.adjust(pvalue))
}



#Figure4B
if(T){
  #FDR
  ttest_result$name <- rownames(ttest_result)
  ttest_result$dis <- (ttest_result$mean_High - ttest_result$mean_Low)
  ttest_result$H_L <- ifelse(ttest_result$dis > 0,"H","L")
  ttest_result$`q-value` <- NA
  ttest_result[which(ttest_result$qvalue > 0.1),'q-value'] <- 'NA'
  ttest_result[which(ttest_result$qvalue < 0.1 & ttest_result$qvalue > 0.05),'q-value'] <- '*'
  ttest_result[which(ttest_result$qvalue < 0.05 & ttest_result$qvalue > 0.01),'q-value'] <- '**'
  ttest_result[which(ttest_result$qvalue < 0.01),'q-value'] <- '***'
  ttest_result$`q-value` <- factor(ttest_result$'q-value',levels = c('NA','*','**','***'))
  
  ttest_result %<>%
    arrange(dis)
  
  plot4B <- ggdotchart(ttest_result, x = "name", y = "dis",
    color = 'q-value', 
    palette = c("#8e9eab", "#ffba08","#e85d04","#d00000"), 
    size = 6, legend.title = 'T-test FDR',
    label = 'H_L', 
    font.label = list(color = "white", size = 9, vjust = 0.5), 
    add = "segments", 
    add.params = list(color = "lightgray",size = 1.5),
    ggtheme = theme_pubr(),
    xlab = '', ylab = '',
    sorting = "descending", 
    rotate = TRUE) +
    scale_y_continuous(breaks = c(-0.05,0,0.05), labels = c('Low_Score','|','High_Score'))+
    theme(legend.position = c(0.90,0.18), axis.ticks.x = element_blank())
}
plot4B


#Figure4C
if(T){
    library(ggraph)
    library(igraph)
    cor_graph <- cor_result %>% 
      mutate(type = ifelse(cor > 0, "pos", "neg"),
        from = "ES") %>% 
      mutate(value = abs(cor)) %>% 
      rownames_to_column(var = "to") %>% 
      select("to","value","pvalue","type")
    cor_graph$type[which(cor_graph$pvalue > 0.05)] <- NA
    
    cor_edges <- readxl::read_excel("./input/IMMUNE_cell_type.xlsx",
      sheet = "Sheet3") %>% 
      left_join(cor_graph)
    
    cor_vertices  <-  data.frame(
      name = c("ES",as.character(cor_edges$to)) , 
      value = c(NA,as.numeric(cor_edges$value)),
      pvalue = c(NA,as.numeric(cor_edges$pvalue)),
      type = c(NA,cor_edges$type)) 
    cor_vertices$type[is.na(cor_vertices$type)] <- "None"
    
    mygraph <- graph_from_data_frame(cor_edges, vertices = cor_vertices, directed = TRUE)
    
    plot4C <- ggraph(mygraph,layout = 'dendrogram', circular = TRUE) + 
      geom_edge_diagonal(aes(colour=..index..)) +
      scale_edge_colour_distiller(palette = "RdPu") +
      geom_node_text(aes(x = x*1.25, y=y*1.25, label=name,color=type),
        size=3, alpha=1) +
      geom_node_point(aes(x = x*1.07, y= y*1.07,fill = type,filter = leaf,size = value),
        shape=21,stroke=0,color = "red",alpha=1) +
      scale_colour_manual(values= c('#0000ff',"#000000",'#ff0000')) +
      scale_fill_manual(values= c('#0000ff',"#000000",'#ff0000')) +
      expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
      theme_void()
}
plot4C


# 10. Figure4D ------------------------------------------------------------


rm(list = setdiff(ls(),c('fpkm_rf_gene','fpkm_survival','tcga_pscore')))
load("./input/rdata/RCircos.rdata")
RCircos.Gene.Label.Data %<>% 
  rename("hgnc_symbol" = "Gene") %>% 
  left_join(name_all) %>% 
  filter(gene_biotype == "protein_coding") %>% 
  rename("Gene" = "hgnc_symbol") %>% 
  select(-"gene_biotype")
library(RCircos)
#Figure4D
tracks.inside <- 6
tracks.outside <- 1
RCircos.Set.Core.Components(cyto.info, chr.exclude=c("chrX", "chrY"),tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area()     
RCircos.Chromosome.Ideogram.Plot() 
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col = 4,track.num = 1, side = "out");
RCircos.Gene.Connector.Plot(RCircos.GeneBand.Label.Data, track.num = 1, side = "in");
RCircos.Gene.Name.Plot(RCircos.GeneBand.Label.Data, name.col = 4,track.num = 2, side = "in");
RCircos.Scatter.Plot(data, data.col=10, track.num=5, side="in",by.fold=1)



# 11. Figure5A-G ----------------------------------------------------------


#GDSC file: 
#Cell_line_RMA_proc_basalExp.txt.gz
#PANCANCER_IC_Tue May 31 11_37_12 2022.csv.gz
#PANCANCER_IC_Tue May 31 07_45_17 2022.csv.gz
#download from https://www.cancerrxgene.org/downloads/bulk_download

#tidying up and get rdata

rm(list = setdiff(ls(),c('fpkm_rf_gene','fpkm_survival','tcga_pscore')))

load("./input/rdata/result_GDSC.rdata")

#Figure5A-G

#set the figure DIR
savepicturedir <- './figure/'
sigDrug_GDSC <- result_GDSC %>% 
  filter((Drug_name %in% c('Doxorubicin','Docetaxel','Paclitaxel',
    'Vincristine','Cisplatin','Fulvestrant','Tamoxifen')))

map(sigDrug_GDSC$Drug_name,function(a){
  GDSC_temp <- drug_GDSC %>% 
    filter(Drug.name == a) %>% 
    left_join(clust_GDSC) %>% 
    select("IC50","score") %>% 
    mutate(IC50 = as.numeric(IC50),
      score = as.numeric(score))
  cor.test <- cor.test(as.numeric(GDSC_temp$IC50),as.numeric(GDSC_temp$score),method = "pearson")
  p = cor.test[["p.value"]]
  estimate = cor.test[["estimate"]][["cor"]]
  theme_set(ggpubr::theme_pubr()+
      theme(legend.position = "top"))
  pic <- ggplot(GDSC_temp,aes(x = score,y = IC50))+
    geom_point()+
    stat_smooth(method = lm,color = "black", fill = "lightgray")+
    labs(x = "ES", y = str_c(a," sensitivity (IC50)"))+
    stat_cor(method = "pearson")
  ggsave(str_c(savepicturedir,a,'.pdf',sep = ""),width = 8, height = 8)
})


# 12.Figure5H -------------------------------------------------------------

rm(list = setdiff(ls(),c('fpkm_rf_gene')))

GSE145325_pdata1 <- read.delim("./input/GSE145325/GSE145325_series_matrix.txt", 
  header= T)[c(1,10),2:59] %>% 
  t() 
colnames(GSE145325_pdata1) <- c("sample_name","group")

GSE145325_pdata1 <- GSE145325_pdata1 %>% 
  as.data.frame() %>% 
  mutate(group = str_split(group,pattern = "[ ]",simplify = T)[,2]) %>% 
  rownames_to_column(var = "experiment_alias") %>% 
  mutate(experiment_alias = str_replace(experiment_alias,pattern = "[.]",replacement = "-"))

GSE145325_pdata <- read.delim("./input/GSE145325/filereport_read_run_PRJNA605185_tsv.txt") %>% 
  left_join(GSE145325_pdata1)

#count data
GSE145325_counts <- read.delim("./input/GSE145325/counts.txt.gz", comment.char="#") %>% 
  filter(!str_detect(Chr,pattern = "chrY")) %>% 
  filter(!str_detect(Chr,pattern = "chrM")) %>% 
  select(-c(2:6)) %>% 
  column_to_rownames(var = "Geneid") 

colnames(GSE145325_counts) <- str_split(str_split(colnames(GSE145325_counts),
  pattern = "[.]",simplify = T)[,7],
  pattern = "_",simplify = T)[,1]


#counts to fpkm

gene_length <- read.delim("./input/GSE145325/counts.txt.gz", comment.char="#")[,c(1,6)] %>% 
  filter(Geneid %in% rownames(GSE145325_counts))

GSE145325_fpkm <- sapply(GSE145325_counts,function(x){
  x = (x / sum(x)) * 1e6
  x = x / gene_length$Length
}) %>% 
  as.data.frame()
rownames(GSE145325_fpkm) <- rownames(GSE145325_counts)


GSE145325_score <- GSE145325_fpkm[fpkm_rf_gene[,1],] %>% 
  t() %>% 
  apply(.,1,function(x){
    tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(fpkm_rf_gene[,2]))
    return(tt)
  }) %>% 
  as.data.frame() %>% 
  rownames_to_column("run_accession") %>% 
  rename(score = ".") %>% 
  mutate(score = score * 100)

GSE145325 <- GSE145325_pdata %>% 
  left_join(GSE145325_score)

#Figure5H
plot5H <- ggplot(GSE145325 , aes(x = group, y = score, fill = group))+
  geom_boxplot(aes(x = group,y = score, fill = group),
    outlier.shape = NA,
    width = .07,
    color = "black")+
  stat_compare_means(label.y = 2.5,label.x = 2,size = 5)+
  geom_point(aes(x = as.numeric(as.factor(group))-0.1,
    y = score,color = group),
    position = position_jitter(width = .05),size = .25, shape = 20) +
  geom_half_violin(aes(fill = group),
    position = position_nudge(x = .15, y = 0),
    width = .2,
    adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
  theme_minimal()+
  scale_color_jco() +
  scale_fill_jco() +
  theme(panel.background=element_rect(colour="black",fill=NA),
    panel.grid=element_blank(),
    legend.position="none",
    legend.background=element_rect(colour=NA,fill=NA))
plot5H
