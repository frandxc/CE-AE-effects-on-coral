library(RColorBrewer)
library(ggpubr)
library(purrr)
library(genefilter)
library(DESeq2)
library(amap)
library(limma)
library(cluster)
library(ggfortify)
library("RColorBrewer")
library("amap")
library(ggplot2)
library("BiocParallel")
library(pheatmap) 
library(dplyr)
library(patchwork)
require(cowplot)
library(gridExtra)
library("Biostrings")
library(patchwork)
library(ggvenn)

#library(greenred)
# BiocManager::install(c("greenred"))
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
##annotation and raw count
datahgraw1 <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/All_Unigene.advance.annotation.csv', header = T, row.names=1) 
datahgraw <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/all_g_count.csv', header = T, row.names=1)
# sym blast results----
file_path <- "D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/symbiont_result/sym"
files <- list.files(path = file_path , full.names = TRUE)
data_sym <- read.csv(paste0(file_path, "/clade_A1.blast.csv"), sep = "\t", header = T)
# 循环读取剩余两个CSV文件并进行rbind拼接
for (i in 2:12) {
  file <- read.csv(files[i], sep = "\t", header = T)
  data_sym <- rbind(data_sym, file)
}
data_sym_filter <- data_sym[data_sym$pident>90,]
# annotation----
datahgraw_NR <- datahgraw1
colnames(datahgraw1)
datahgraw_NR$Nr.annotation= gsub("PREDICTED: ", "", x =c(datahgraw_NR$Nr.annotation))
datahgraw_NR$Nr.annotation= gsub("LOW QUALITY PROTEIN: ", "", x =c(datahgraw_NR$Nr.annotation))
datahgraw_NR$Nr.annotation= gsub("\\[.*]","", x =c(datahgraw_NR$Nr.annotation))#delete [Acropora digitifera] and [.*]
datahgraw_NR$Nr.annotation <- gsub(",.*", "", datahgraw_NR$Nr.annotation)
conditionshg <- factor(c(rep("Control", 4), rep("Ca(OH)2", 4), rep("CO2", 4), rep("CO2+Ca(OH)2", 4), rep("CO2+NaOH", 4)))

# grepl: coal gene filter
datahgraw_r <- datahgraw1
datahg1<- merge(datahgraw[,c(1:4,17:20,5:16)],datahgraw_r,by="row.names")
rownames(datahg1)<- datahg1[,1]
datahg1<- datahg1[,-1]
datahg<- datahg1[,1:20]
datahg[,5]<- rowMeans(datahg1[,6:8])#as the 5th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg[,18]<- rowMeans(datahg1[,19:20])#as the 18th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg <-round(as.matrix(datahg))
head(datahg)
# annotation exctraction
h_annot <- datahgraw1[,1:37]

#DEG----------  
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
samplehg_rowname <- colnames(datahg[,1:20])
samplehg <- data.frame(conditionshg)
rownames(samplehg) <- samplehg_rowname
###matix
ddsFullCountTablehg <- DESeqDataSetFromMatrix(countData = datahg,colData = samplehg, design= ~ conditionshg)
###filter
keep <- rowMeans(counts(ddsFullCountTablehg)) >= 10 # 对counts数进行初步的过滤
dds_hg <- ddsFullCountTablehg[keep,]
#####deseq
ddshg <- DESeq(dds_hg)

###gene normalization data---------
rldhg <- rlogTransformation(ddshg)  ## 得到经过DESeq2软件normlization 
normalized_all=data.frame(assay(rldhg))
filtered_coral <- subset(datahgraw1, grepl("\\[Acropora digitifera\\]|\\[Exaiptasia pallida\\]|\\[Acropora millepora\\] |\\[Nematostella vectensis\\]", Nr.annotation))
normalized_hg <- data.frame(normalized_all[rownames(normalized_all) %in% rownames(filtered_coral),])
normalized_hg_coral <- normalized_hg 
normalized_hg_sym <- data.frame(normalized_all[rownames(normalized_all) %in% data_sym_filter$qaccver,]) 
normalized_all= rbind(normalized_hg_coral,normalized_hg_sym)
normalized_all <- rbind(normalized_hg_coral, normalized_hg_sym)
normalized_all$group <- rep(c("coral", "sym"), times = c(nrow(normalized_hg_coral), nrow(normalized_hg_sym))) # Create a new column for grouping

#pheatmap(normalized_hg,scale='row',border_color=NA,cluster_cols = FALSE,cluster_rows = F,legend =FALSE, labels_row=c(""))

boxplot(normalized_hg, col = rainbow(ncol(normalized_hg)*1.2),main="expression value",las=2)

## negative and positive DEGs-----------
hgConCaOH <- results(ddshg, contrast = c("conditionshg","Ca(OH)2","Control"))
hgConCO2 <- results(ddshg, contrast = c("conditionshg","CO2","Control"))
hgConCO2CaOH <- results(ddshg, contrast = c("conditionshg","CO2+Ca(OH)2","Control"))
hgConCO2NaOH <- results(ddshg, contrast = c("conditionshg","CO2+NaOH","Control"))

ConCaOHlist <- data.frame(rowname=hgConCaOH@rownames, hgConCaOH@listData)
rownames(ConCaOHlist) <-ConCaOHlist$rowname
mcaoh =merge(normalized_hg,ConCaOHlist , by="row.names")
# PLOT: gene UP AND DOWN  NUMBER ------
up_CaOH=nrow(subset(hgConCaOH,log2FoldChange>2&padj<0.01))
dwon_CaOH=nrow(subset(hgConCaOH,log2FoldChange< -2 & padj<0.01))
up_CO2=nrow(subset(hgConCO2,log2FoldChange>2&padj<0.01))
dwon_CO2=nrow(subset(hgConCO2,log2FoldChange< -2 & padj<0.01))
up_CO2CaOH=nrow(subset(hgConCO2CaOH,log2FoldChange>2&padj<0.01))
dwon_CO2CaOH=nrow(subset(hgConCO2CaOH,log2FoldChange< -2 & padj<0.01))
up_CO2NaOH=nrow(subset(hgConCO2NaOH ,log2FoldChange>2&padj<0.01))
dwon_CO2NaOH=nrow(subset(hgConCO2NaOH ,log2FoldChange< -2 & padj<0.01))
dat1 <- data.frame(up_CaOH,dwon_CaOH,up_CO2,dwon_CO2,up_CO2CaOH,dwon_CO2CaOH,up_CO2NaOH,dwon_CO2NaOH) 
tdat1 <- data.frame(treat=colnames(dat1),number=t(dat1))
tdat1$treat <- factor(tdat1$treat, level= c(tdat1$treat),ordered=TRUE)
pdf("DEG_COAL_number.pdf", width=3, height=3)
ggplot(tdat1,aes(treat,number))+geom_bar(stat="identity",fill=c("blue","red","blue","red","blue","red","blue","red")) +theme_classic()+  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
# upregulate mRNA-----
upmRNA_CO2<-rownames(subset(hgConCO2,log2FoldChange>2&padj<0.01))  #HC upregulate  
upmRNA_CO2CaOH<-rownames(subset(hgConCO2CaOH,log2FoldChange>2&padj<0.01))  #HC upregulate  CO2CaOH2_up$name  #HCHA1 upregulate
upmRNA_CO2NaOH<-rownames(subset(hgConCO2NaOH,log2FoldChange>2&padj<0.01))  #HC upregulate  CO2NaOH_up$name  #HCHA2 upregulate
upmRNA_CaOH<-rownames(subset(hgConCaOH,log2FoldChange>2&padj<0.01))  #HC upregulate  CO2NaOH_up$name  #HCHA2 upregulate

upmRNA_reverse<-list(`CO2_mRNA`=upmRNA_CO2, 
                   `CO2CaOH_mRNA`= upmRNA_CO2CaOH,
                   `CO2NaOH_mRNA`= upmRNA_CO2NaOH)
#venn plot----
ggvenn_upRNA<- ggvenn(upmRNA_reverse,
       show_percentage =F,
       stroke_color = "NA",fill_color =  c("#ffb2b2", "#b2e7cb", "#ffcc99"),
       set_name_color = c("#ff0000", "#4a9b83", "#0000ff"))
# downregulate mRNA------
downmRNA_CO2<-rownames(subset(hgConCO2,log2FoldChange< -2&padj<0.01))  #HC downregulate  
downmRNA_CO2CaOH<-rownames(subset(hgConCO2CaOH,log2FoldChange< -2&padj<0.01))  #HC downregulate  CO2CaOH2_down$name  #HCHA1 downregulate
downmRNA_CO2NaOH<-rownames(subset(hgConCO2NaOH,log2FoldChange< -2&padj<0.01))  #HC downregulate  CO2NaOH_down$name  #HCHA2 downregulate
downmRNA_CaOH<-rownames(subset(hgConCaOH,log2FoldChange< -2&padj<0.01))  #HC downregulate  CO2CaOH2_down$name  #HCHA1 downregulate
downmRNA_reverse<-list(`CO2_mRNA`=downmRNA_CO2, 
                     `CO2CaOH_mRNA`= downmRNA_CO2CaOH,
                     `CO2NaOH_mRNA`= downmRNA_CO2NaOH)
#venn plot----
ggvenn_downRNA<- ggvenn(downmRNA_reverse,
       show_percentage =F,
       stroke_color = "NA",fill_color =  c("#ffb2b2", "#b2e7cb", "#ffcc99"),
       set_name_color = c("#ff0000", "#4a9b83", "#0000ff"))
library(ggplot2)
theme_set(theme_classic())

# PCA plot------

pca_result <- prcomp(t(data.frame(normalized_hg)), scale. = TRUE) #normalized_hg is normalized host mRNA
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Category = conditionshg)

# 计算PC1和PC2的百分比
pca_var_ratio <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

theme_set(theme_classic())
PCA_DEG<- ggplot(pca_data, aes(x = PC1, y = PC2, color = Category)) +
  geom_point(size = 3, alpha = 0.6) +
  labs(x = paste("PC1 (", round(pca_var_ratio[1], 2), "%)"), 
       y = paste("PC2 (", round(pca_var_ratio[2], 2), "%)"),
       title = "PCA Plot") +
  scale_color_manual(values = c('red', 'blue', 'green', 'purple', 'brown'), labels = unique(conditionshg)) +
  theme(legend.title = element_blank())

ggsave("PCA_DEG.pdf",PCA_DEG,w=4,h=4)


# DEP-------------------
protein_data_r <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/protein_LFQ.csv', header = T) 

datahgraw <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/All_Unigene.advance.annotation.csv', header = T) 
datahgraw_n<- data.frame(datahgraw[,1])
colnames(datahgraw_n) <- c("ID")
h_annot <- datahgraw[,1:37]
##merge with host gene to find host protein
protein_data_coral <- protein_data_r 
  
protein_data<-merge(protein_data_coral,datahgraw_n, by.x='Protein',by.y='ID')
# Generate experimental design
library("stringr")
library("DEP")
library("tidyr")
library(dplyr)
library("SummarizedExperiment")
replicates=4
experimental_design <- tibble(
  label =str_split_fixed(colnames(protein_data)[grep("LFQ.", colnames(protein_data))],'y.',n=2)[,2],
  condition = conditionshg,
  replicate = rep(1:replicates, 5))
# Generate a SummarizedExperiment object using an experimental design
protein_data$Protein %>% duplicated() %>% any()
# protein_data %>% group_by(Protein) %>% summarize(frequency = n()) %>% 
#   arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(protein_data, "Protein", "FastaHeaders", delim = ";")
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
data_filt <- filter_missval(data_se, thr = 0)
data_filt2 <- filter_missval(data_se, thr = 1)
# Normalize the data
data_norm <- normalize_vsn(data_filt)
df1=data.frame(data_norm@assays@data@listData)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)
impute(data_norm, fun = "MinProb")
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_imp_norm=data.frame(data_imp@assays@data@listData)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)
####normalized data-----
data_imp_all <-data.frame(data_imp@assays@data)
data_imp_all$Protein <- rownames(data_imp_all)
data_imp_all_m <- merge(data_imp_all, protein_data_r[,c(1,34:50)], by = "Protein")
data_imp_coral <- subset(data_imp_all_m, grepl("\\[Acropora digitifera\\]|\\[Exaiptasia pallida\\]|\\[Acropora millepora\\] |\\[Nematostella vectensis\\]\\[Pocillopora damicornis\\]", Nr.annotation))
data_imp_sym <- data_imp_all_m[data_imp_all_m$Protein %in% data_sym_filter$qaccver,]
data_imp_all1 <- rbind(data_imp_coral,data_imp_sym) [,-3]
data_imp_all1$group <- rep(c("coral", "sym"), times = c(nrow(data_imp_coral), nrow(data_imp_sym)))

# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Control")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.2))
dep_all  <- data.frame(dep@assays@data@listData[[1]])
dep_coral  <- data.frame(dep_all[rownames(dep_all) %in% data_imp_coral$Protein,])

plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

#PCA plot------
pca_DAP <- prcomp(t(data.frame(data_imp_coral[,4:23])), scale. = TRUE)
pca_data_DAP <- data.frame(PC1 = pca_DAP$x[,1], PC2 = pca_DAP$x[,2], Category = conditionshg)

# 计算PC1和PC2的百分比
pca_DAP_ratio <- pca_DAP$sdev^2 / sum(pca_DAP$sdev^2) * 100

theme_set(theme_classic())
PCA_DAP<- ggplot(pca_data_DAP, aes(x = PC1, y = PC2, color = Category)) +
  geom_point(size = 3, alpha = 0.6) +
  labs(x = paste("PC1 (", round(pca_DAP_ratio[1], 2), "%)"), 
       y = paste("PC2 (", round(pca_DAP_ratio[2], 2), "%)"),
       title = "PCA_DAP Plot") +
  scale_color_manual(values = c('red', 'blue', 'green', 'purple', 'brown'), labels = unique(conditionshg)) +
  theme(legend.title = element_blank())

ggsave("PCA_DAP.pdf",PCA_DAP,w=4,h=4)
#all plot -----
library(cowplot)
# 组合所有图形并设置主图的尺寸

plot_grid(PCA_DEG, PCA_DAP,
          ncol = 2, align = "v",
          labels = c("a", "b",
                     label_size = 14,
                     rel_heights = c(1, 1),
                     rel_widths = c(1, 1)))
ggsave("PCA_DEG_DAP.pdf",w=8,h=4)


plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
data_results <- get_results(dep)
data_results_coral  <- data_results[data_results$name %in% data_imp_coral$Protein,]
# Number of significant proteins
data_results_coral %>% filter(significant) %>% nrow()
# Column names of the results table
colnames(data_results)


##positive and negative DEPs----

CaOH2_up <- data_results[data_results$Ca.OH.2_vs_Control_p.adj<0.05&data_results$Ca.OH.2_vs_Control_ratio >=1.2,][grep("^Ca.OH.2_vs|name", colnames(data_results))]
CaOH2_down <- data_results[data_results$Ca.OH.2_vs_Control_p.adj<0.05&data_results$Ca.OH.2_vs_Control_ratio<=-1.2,][grep("^Ca.OH.2_vs|name", colnames(data_results))]

CO2_up <- data_results[data_results$CO2_vs_Control_p.adj<0.05&data_results$CO2_vs_Control_ratio>=1.2,][grep("^CO2_vs|name", colnames(data_results))]
CO2_down <- data_results[data_results$CO2_vs_Control_p.adj<0.05&data_results$CO2_vs_Control_ratio<=-1.2,][grep("^CO2_vs|name", colnames(data_results))]

CO2CaOH2_up <- data_results[data_results$CO2.Ca.OH.2_vs_Control_p.adj<0.05&data_results$CO2.Ca.OH.2_vs_Control_ratio>=1.2,][grep("^CO2.Ca.OH.2_vs|name", colnames(data_results))]
CO2CaOH2_down <- data_results[data_results$CO2.Ca.OH.2_vs_Control_p.adj<0.05&data_results$CO2.Ca.OH.2_vs_Control_ratio<=-1.2,][grep("^CO2.Ca.OH.2_vs|name", colnames(data_results))]

CO2NaOH_up <- data_results[data_results$CO2.NaOH_vs_Control_p.adj<0.05&data_results$CO2.NaOH_vs_Control_ratio>=1.2,][grep("^CO2.NaOH_vs|name", colnames(data_results))]
CO2NaOH_down <- data_results[data_results$CO2.NaOH_vs_Control_p.adj <0.05&data_results$CO2.NaOH_vs_Control_ratio<=-1.2,][grep("^CO2.NaOH_vs|name", colnames(data_results))]

datPr <- data.frame(treat=c(tdat1$treat),V1=c(nrow(CaOH2_up),nrow(CaOH2_down),nrow(CO2_up ),nrow(CO2_down),nrow(CO2CaOH2_up),nrow(CO2CaOH2_down),nrow(CO2NaOH_up),nrow(CO2NaOH_down)))
datPr$treat<- factor(tdat1$treat, level= c(tdat1$treat),ordered=TRUE)
pdf("DEP_COAL_number.pdf", width=3, height=3)
ggplot(datPr,aes(treat,V1))+geom_bar(stat="identity",fill=c("blue","red","blue","red","blue","red","blue","red")) +theme_classic()+  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

ANNO_CO2 <- datahgraw1[rownames(datahgraw1) %in% CaOH2_up$name,]
#Fig 2A VENN plot of mRNA and protein ---------------------

library(ggvenn)
library(VennDiagram)
genelist_CaOH= data.frame(rownames(hgConCaOH),hgConCaOH$log2FoldChange,hgConCaOH$pvalue)
names(genelist_CaOH)=c('GID','FC','pvalue')
genelist_CaOH= subset(genelist_CaOH,abs(genelist_CaOH$FC)>2&hgConCaOH$pvalue<0.05)
#mRNA caOH## positive 
genelist_CaOH_pr= data_results[data_results$Ca.OH.2_vs_Control_p.val<0.05&abs(data_results$Ca.OH.2_vs_Control_ratio)>1.2,][grep("^Ca.OH.2_vs|name", colnames(data_results))]
genelist_CaOH_p=data.frame(GID=genelist_CaOH_pr$name,FC=genelist_CaOH_pr$Ca.OH.2_vs_Control_ratio,pvalue=genelist_CaOH_pr$Ca.OH.2_vs_Control_p.val)

names(genelist_CaOH)=c('GID','FC','pvalue')
names(genelist_CaOH_p)=c('GID','FC','pvalue')

CaOH_VENN<-list(`mRNA`=genelist_CaOH$GID, 
                `protein`= genelist_CaOH_p$GID)


#mRNA CO2## positive 
genelist_CO2= data.frame(rownames(hgConCO2),hgConCO2$log2FoldChange,hgConCO2$pvalue)
names(genelist_CO2)=c('GID','FC','pvalue')
genelist_CO2= subset(genelist_CO2,abs(genelist_CO2$FC)>2&hgConCO2$pvalue<0.05)
#mRNA CO2## positive 
genelist_CO2_pr= data_results[data_results$CO2_vs_Control_p.val<0.05&abs(data_results$CO2_vs_Control_ratio)>1.2,][grep("^CO2_vs|name", colnames(data_results))]
genelist_CO2_p=data.frame(GID=genelist_CO2_pr$name,FC=genelist_CO2_pr$CO2_vs_Control_ratio,pvalue=genelist_CO2_pr$CO2_vs_Control_p.val)

names(genelist_CO2)=c('GID','FC','pvalue')


CO2_VENN<-list(`mRNA`=genelist_CO2$GID, 
               `protein`= genelist_CO2_p$GID)
#mRNA CO2CaOH2## positive 
genelist_CO2CaOH2= data.frame(rownames(hgConCO2CaOH),hgConCO2CaOH$log2FoldChange,hgConCO2CaOH$pvalue)
names(genelist_CO2CaOH2)=c('GID','FC','pvalue')
genelist_CO2CaOH2= subset(genelist_CO2CaOH2,abs(genelist_CO2CaOH2$FC)>2&hgConCO2CaOH$pvalue<0.05)
#mRNA CO2CaOH2## positive 
genelist_CO2CaOH2_pr= data_results[data_results$CO2.Ca.OH.2_vs_Control_p.val<0.05&abs(data_results$CO2.Ca.OH.2_vs_Control_ratio)>1.2,][grep("^CO2.Ca.OH.2_vs|name", colnames(data_results))]
genelist_CO2CaOH2_p=data.frame(GID=genelist_CO2CaOH2_pr$name,FC=genelist_CO2CaOH2_pr$CO2.Ca.OH.2_vs_Control_ratio,pvalue=genelist_CO2CaOH2_pr$CO2.Ca.OH.2_vs_Control_p.val)

names(genelist_CO2CaOH2)=c('GID','FC','pvalue')


CO2CaOH2_VENN<-list(`mRNA`=genelist_CO2CaOH2$GID, 
                    `protein`= genelist_CO2CaOH2_p$GID)

#mRNA CO2NaOH## positive 
genelist_CO2NaOH= data.frame(rownames(hgConCO2NaOH),hgConCO2NaOH$log2FoldChange,hgConCO2NaOH$pvalue)
names(genelist_CO2NaOH)=c('GID','FC','pvalue')
genelist_CO2NaOH= subset(genelist_CO2NaOH,abs(genelist_CO2NaOH$FC)>2&hgConCO2NaOH$pvalue<0.05)
#mRNA CO2NaOH## positive 
genelist_CO2NaOH_pr= data_results[data_results$CO2.NaOH_vs_Control_p.val<0.05&abs(data_results$CO2.NaOH_vs_Control_ratio)>1.2,][grep("^CO2.NaOH_vs|name", colnames(data_results))]
genelist_CO2NaOH_p=data.frame(GID=genelist_CO2NaOH_pr$name,FC=genelist_CO2NaOH_pr$CO2.NaOH_vs_Control_ratio,pvalue=genelist_CO2NaOH_pr$CO2.NaOH_vs_Control_p.val)

names(genelist_CO2NaOH)=c('GID','FC','pvalue')


CO2NaOH_VENN<-list(`mRNA`=genelist_CO2NaOH$GID, 
                   `protein`= genelist_CO2NaOH_p$GID)
#
plot_list = list() 


plot_list[[1]]=ggvenn(CaOH_VENN,
                      show_percentage =F,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))

plot_list[[2]]=ggvenn(CO2_VENN,
                      show_percentage =F,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))

plot_list[[3]]=ggvenn(CO2CaOH2_VENN,
                      show_percentage =F,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))
plot_list[[4]]=ggvenn(CO2NaOH_VENN,
                      show_percentage =F,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))

pdf("Fig 2A gene_protein_venn.pdf",width=20,height=30) 
pushViewport(viewport(layout = grid.layout(2,2)))
print(plot_list[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col =1))
print(plot_list[[2]], vp = viewport(layout.pos.row = 2, layout.pos.col =1))
print(plot_list[[3]], vp = viewport(layout.pos.row = 1, layout.pos.col =2))
print(plot_list[[4]], vp = viewport(layout.pos.row = 2, layout.pos.col =2))
dev.off() 


# upregulate protein------
genelist_CO2p<-CO2_up$name  #HC upregulate  
genelist_CO2CaOHp<-CO2CaOH2_up$name  #HCHA1 upregulate
genelist_CO2NaOHp<-CO2NaOH_up$name  #HCHA2 upregulate
genelist_CaOHp<-CaOH2_up$name  #HCHA1 upregulate
up_p_reverse<-list(`CO2_P`=CO2_up$name, 
                   `CO2CaOH_P`= CO2CaOH2_up$name,
                   `CO2NaOH_P`= CO2NaOH_up$name)
#venn plot----
ggvenn_up_p<- ggvenn(up_p_reverse,
       show_percentage =F,
       stroke_color = "NA",fill_color =  c("#ffb2b2", "#b2e7cb", "#ffcc99"),
       set_name_color = c("#ff0000", "#4a9b83", "#0000ff"))


#downregulate mRNA
genelist_CO2d<-CO2_down$name
genelist_CO2CaOHd<-CO2CaOH2_down$name
genelist_CO2NaOHd<-CO2NaOH_down$name
gene_reverse_d<-list(`CO2`=genelist_CO2d, 
                     `CO2CaOH`= genelist_CO2CaOHd,
                     `CO2NaOH`= genelist_CO2NaOHd)
inter_down <- get.venn.partitions(gene_reverse_d)



#upregulated HC heatmap-----
inter_up <- get.venn.partitions(up_p_reverse)
up_HC <- unlist(inter_up[7,5])
p_up_HC_heatmap<-data_imp_coral[data_imp_coral$Protein %in% up_HC, ][,c(1,4:23)]
p_up_HC_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_up_HC_heatmap$Protein,]

# upregulated HCHA heatmap-----

up_HCHA <- unlist(inter_up[2,5])
p_up_HCHA_heatmap<-data_imp_coral[data_imp_coral$Protein %in% up_HCHA, ][,c(1,4:23)]
p_up_HCHA_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_up_HCHA_heatmap$Protein,]

# upregulated HC_HCHA------
up_HC_HCHA<- unlist(inter_up[1,5])
p_up_HC_HCHA_heatmap<-data_imp_coral[data_imp_coral$Protein %in% up_HC_HCHA, ][,c(1,4:23)]
p_up_HC_HCHA_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_up_HC_HCHA_heatmap$Protein,]

# downregulated protein-----
genelist_CO2p<-CO2_down$name  #HC downregulate  
genelist_CO2CaOHp<-CO2CaOH2_down$name  #HCHA1 downregulate
genelist_CO2NaOHp<-CO2NaOH_down$name  #HCHA2 downregulate
down_p_reverse<-list(`CO2_P`=CO2_down$name, 
                   `CO2CaOH_P`= CO2CaOH2_down$name,
                   `CO2NaOH_P`= CO2NaOH_down$name)
#venn plot----
ggvenn_down_p<- ggvenn(down_p_reverse,
       show_percentage =F,
       stroke_color = "NA",fill_color =  c("#ffb2b2", "#b2e7cb", "#ffcc99"),
       set_name_color = c("#ff0000", "#4a9b83", "#0000ff"))

# all venn upset plot----
# 
# (ggvenn_upRNA+ggvenn_up_p)/(ggvenn_downRNA+ggvenn_down_p)
# upmRNA_p_reverse<- append(upmRNA_reverse, up_p_reverse)  
# downmRNA_p_reverse<- append(downmRNA_reverse, down_p_reverse)  
# 
# upmRNA_p_reverse<- append(upmRNA_reverse, up_p_reverse)  

library(UpSetR)

up_mrna_p <- list(up_CaOH_mRNA=rownames(subset(hgConCaOH,log2FoldChange>2&padj<0.01)),
                  up_CO2_mRNA=rownames(subset(hgConCO2,log2FoldChange>2&padj<0.01)),
                  up_CO2CaOH_mRNA=rownames(subset(hgConCO2CaOH,log2FoldChange>2&padj<0.01)),  #HC upregulate  CO2CaOH2_up$name  #HCHA1 upregulate
                  up_CO2NaOH_mRNA=rownames(subset(hgConCO2NaOH,log2FoldChange>2&padj<0.01)),
                  up_CaOH_p=CaOH2_up$name,  #HCHA1 upregulate
                  up_CO2_p=CO2_up$name,
                  up_CO2CaOH_p=CO2CaOH2_up$name,  #HCHA2 upregulate
                  up_CO2NaOH_p=CO2NaOH_up$name)  #HCHA1 upregulate
down_mrna_p <- list(down_CaOH_mRNA=rownames(subset(hgConCaOH,log2FoldChange>2&padj<0.01)),
                    down_CO2_mRNA=rownames(subset(hgConCO2,log2FoldChange>2&padj<0.01)),
                  down_CO2CaOH_mRNA=rownames(subset(hgConCO2CaOH,log2FoldChange>2&padj<0.01)),  #HC downregulate  CO2CaOH2_down$name  #HCHA1 downregulate
                  down_CO2NaOH_mRNA=rownames(subset(hgConCO2NaOH,log2FoldChange>2&padj<0.01)),
                  down_CaOH_p=CaOH2_down$name,  #HCHA1 downregulate
                  down_CO2_p=CO2_down$name,
                  down_CO2CaOH_p=CO2CaOH2_down$name,  #HCHA2 downregulate
                  down_CO2NaOH_p=CO2NaOH_down$name)  #HCHA1 downregulate
upsetData=fromList(up_mrna_p)

setbarcolor <- c("#2e409a", "#942d8d", "#d75427", "#006b7b", "#4da0a0", "#9b3a74", "#FF5733", "#FFC300")
desired_order1 <- c("up_CO2NaOH_mRNA", "up_CO2CaOH_mRNA","up_CO2_mRNA", "up_CaOH_mRNA",
                    "up_CO2NaOH_p", "up_CO2CaOH_p","up_CO2_p","up_CaOH_p")
upset_up<-upset(upsetData, sets = desired_order1, nset = 8, order.by = c('degree','freq'), decreasing = c(F, T),#排序
      sets.bar.color = setbarcolor, keep.order = TRUE,#让集合按照 sets 参数中指定的出现的顺序排列
      mb.ratio = c(0.4,0.6))
inter_up_mrna_p<- get.venn.partitions(up_mrna_p)

downsetData=fromList(down_mrna_p)
desired_order2 <- c("down_CO2NaOH_mRNA", "down_CO2CaOH_mRNA","down_CO2_mRNA", "down_CaOH_mRNA",
                    "down_CO2NaOH_p", "down_CO2CaOH_p","down_CO2_p","down_CaOH_p")
upset_down<- upset(downsetData, sets = desired_order2, nset = 8, order.by = c('degree','freq'), decreasing = c(F, T),#排序
                   sets.bar.color = setbarcolor, keep.order = TRUE,#让集合按照 sets 参数中指定的出现的顺序排列
                   mb.ratio = c(0.4,0.6))
upset_down1 <- ggplotify::as.ggplot(upset_down)
upset_up1 <- ggplotify::as.ggplot(upset_up)
upset_both<- plot_grid(upset_up1,upset_down1,ncol = 1, align = "hv"  # "hv" stands for horizontal and vertical alignment
)
inter_down_mrna_p<- get.venn.partitions(down_mrna_p)

datahgraw1[rownames(datahgraw1) %in% unlist(inter_down_mrna_p[9,]$..values..),]
datahgraw1[rownames(datahgraw1) %in% unlist(inter_down_mrna_p[17,]$..values..),]
ggsave("upset_both.pdf",upset_both, w=10,h=5)


#downregulation HC-----
inter_down <- get.venn.partitions(gene_reverse_d)
inter_down_HC <- unlist(inter_down[7,5])
p_down_HC_heatmap<-data_imp_coral[data_imp_coral$Protein %in% inter_down_HC, ][,c(1,4:23)]
p_down_HC_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_down_HC_heatmap$Protein,]

#downregulation HC_HCHA -----
inter_down_HC_HCHA<- unlist(inter_down[1,5])
p_down_HC_HCHA_heatmap<-data_imp_coral[data_imp_coral$Protein %in% inter_down_HC_HCHA, ][,c(1,4:23)]
p_down_HC_HCHA_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_down_HC_HCHA_heatmap$Protein,]

#downregulation HCHA-----
inter_down_HCHA <- unlist(inter_down[2,5])
p_down_HCHA_heatmap<-data_imp_coral[data_imp_coral$Protein %in% inter_down_HCHA, ][,c(1,4:23)]
p_down_HCHA_heatmap[,22:(ncol(datahgraw_NR)+22)]<-datahgraw_NR[rownames(datahgraw_NR) %in% p_down_HCHA_heatmap$Protein,]



# plot upregulated-----------
library(Peptides) #calculate pI
library(rentrez) #extract sequence

#HC UP 
HC_sequences <- entrez_fetch(db = "protein", id = p_up_HC_heatmap$Nr.ID, rettype = "fasta")
write(HC_sequences, file="HC_sequence.fasta", sep="\n")
pI_HC <- readBStringSet("HC_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HC_name1 <- names(pI_HC) %>% as.data.frame()
pI_HC_data1 <- as.data.frame(pI_HC)
pI_HC_data21 <- data.frame(ID= pI_HC_name1$.,seq=pI_HC_data1$x, ID2= substr( pI_HC_name1$., 1, 14))
pI_HC_seq1=data.frame(pI_HC_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HC_data21$seq),1)) 
p_up_HC_heatmap1<- p_up_HC_heatmap[,2:21]
p_up_HC_heatmap1$pI<- pI_HC_seq1$PI1
p_up_HC_heatmap1$ann<- p_up_HC_heatmap$Nr.annotation
p_up_HC_heatmap1$Protein<- p_up_HC_heatmap$Protein
p_up_HC_heatmap2 <- p_up_HC_heatmap1[order(p_up_HC_heatmap1$pI), ]
p_up_HC_heatmap2_mRNA <-  normalized_hg[p_up_HC_heatmap2$Protein,] 
#protein heatmap
up_HC_heat_p=pheatmap(
  p_up_HC_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
up_HC_heat_mRNA=pheatmap(
  p_up_HC_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)

up_HC_heat_p_bar <- ggplot(p_up_HC_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_minimal() +  xlim(0,15)+
geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线
 

#HC and HCHA UP 
HC_HCHA_sequences <- entrez_fetch(db = "protein", id = p_up_HC_HCHA_heatmap$Nr.ID, rettype = "fasta")
write(HC_HCHA_sequences, file="HC_HCHA_sequence.fasta", sep="\n")
pI_HC_HCHA <- readBStringSet("HC_HCHA_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HC_HCHA_name1 <- names(pI_HC_HCHA) %>% as.data.frame()
pI_HC_HCHA_data1 <- as.data.frame(pI_HC_HCHA)
pI_HC_HCHA_data21 <- data.frame(ID= pI_HC_HCHA_name1$.,seq=pI_HC_HCHA_data1$x, ID2= substr( pI_HC_HCHA_name1$., 1, 14))
pI_HC_HCHA_seq1=data.frame(pI_HC_HCHA_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HC_HCHA_data21$seq),1)) 
p_up_HC_HCHA_heatmap1<- p_up_HC_HCHA_heatmap[,2:21]
p_up_HC_HCHA_heatmap1$pI<- pI_HC_HCHA_seq1$PI1
p_up_HC_HCHA_heatmap1$ann<- p_up_HC_HCHA_heatmap$Nr.annotation
p_up_HC_HCHA_heatmap1$Protein<- p_up_HC_HCHA_heatmap$Protein
p_up_HC_HCHA_heatmap2 <- p_up_HC_HCHA_heatmap1[order(p_up_HC_HCHA_heatmap1$pI), ]
p_up_HC_HCHA_heatmap2_mRNA <-  normalized_hg[p_up_HC_HCHA_heatmap2$Protein,] 
#protein heatmap
up_HC_HCHA_heat_p=pheatmap(
  p_up_HC_HCHA_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
up_HC_HCHA_heat_mRNA=pheatmap(
  p_up_HC_HCHA_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
up_HC_HCHA_heat_p_bar <- ggplot(p_up_HC_HCHA_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_minimal() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线


#HCHA UP 
HCHA_sequences <- entrez_fetch(db = "protein", id = p_up_HCHA_heatmap$Nr.ID, rettype = "fasta")
write(HCHA_sequences, file="HCHA_sequence.fasta", sep="\n")
pI_HCHA <- readBStringSet("HCHA_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HCHA_name1 <- names(pI_HCHA) %>% as.data.frame()
pI_HCHA_data1 <- as.data.frame(pI_HCHA)
pI_HCHA_data21 <- data.frame(ID= pI_HCHA_name1$.,seq=pI_HCHA_data1$x, ID2= substr( pI_HCHA_name1$., 1, 14))
pI_HCHA_seq1=data.frame(pI_HCHA_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HCHA_data21$seq),1)) 
p_up_HCHA_heatmap1<- p_up_HCHA_heatmap[,2:21]
p_up_HCHA_heatmap1$pI<- pI_HCHA_seq1$PI1
p_up_HCHA_heatmap1$ann<- p_up_HCHA_heatmap$Nr.annotation
p_up_HCHA_heatmap1$Protein<- p_up_HCHA_heatmap$Protein
p_up_HCHA_heatmap2 <- p_up_HCHA_heatmap1[order(p_up_HCHA_heatmap1$pI), ]
p_up_HCHA_heatmap2_mRNA <-  normalized_hg[p_up_HCHA_heatmap2$Protein,] 
#protein heatmap
up_HCHA_heat_p=pheatmap(
  p_up_HCHA_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
up_HCHA_heat_mRNA=pheatmap(
  p_up_HCHA_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
up_HCHA_heat_p_bar <- ggplot(p_up_HCHA_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_minimal() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线
all_plots <- plot_grid(up_HC_heat_mRNA$gtable, up_HC_heat_p$gtable, 
                        up_HC_HCHA_heat_mRNA$gtable, up_HC_HCHA_heat_p$gtable, 
                        up_HCHA_heat_mRNA$gtable, up_HCHA_heat_p$gtable, ncol = 2, align = "hv"  # "hv" stands for horizontal and vertical alignment
)
bar_plot<-plot_grid(up_HC_heat_p_bar,up_HC_HCHA_heat_p_bar,up_HCHA_heat_p_bar,ncol = 1,align = "hv")
ggsave("pI_UP_combined_heatmaps.pdf", all_plots, width =6, height = 15)
ggsave("pI_UP_combined_bar.pdf", bar_plot, width = 9, height = 15)

# plot downregulation------
HC_sequences <- entrez_fetch(db = "protein", id = p_down_HC_heatmap$Nr.ID, rettype = "fasta")
write(HC_sequences, file="HC_sequence.fasta", sep="\n")
pI_HC <- readBStringSet("HC_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HC_name1 <- names(pI_HC) %>% as.data.frame()
pI_HC_data1 <- as.data.frame(pI_HC)
pI_HC_data21 <- data.frame(ID= pI_HC_name1$.,seq=pI_HC_data1$x, ID2= substr( pI_HC_name1$., 1, 14))
pI_HC_seq1=data.frame(pI_HC_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HC_data21$seq),1)) 
p_down_HC_heatmap1<- p_down_HC_heatmap[,2:21]
p_down_HC_heatmap1$pI<- pI_HC_seq1$PI1
p_down_HC_heatmap1$ann<- p_down_HC_heatmap$Nr.annotation
p_down_HC_heatmap1$Protein<- p_down_HC_heatmap$Protein
p_down_HC_heatmap2 <- p_down_HC_heatmap1[order(p_down_HC_heatmap1$pI), ]
p_down_HC_heatmap2_mRNA <-  normalized_hg[p_down_HC_heatmap2$Protein,] 
#protein heatmap
down_HC_heat_p=pheatmap(
  p_down_HC_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
down_HC_heat_mRNA=pheatmap(
  p_down_HC_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
down_HC_heat_p_bar <- ggplot(p_down_HC_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_minimal() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线


#HC and HCHA down 
HC_HCHA_sequences <- entrez_fetch(db = "protein", id = p_down_HC_HCHA_heatmap$Nr.ID, rettype = "fasta")
write(HC_HCHA_sequences, file="HC_HCHA_sequence.fasta", sep="\n")
pI_HC_HCHA <- readBStringSet("HC_HCHA_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HC_HCHA_name1 <- names(pI_HC_HCHA) %>% as.data.frame()
pI_HC_HCHA_data1 <- as.data.frame(pI_HC_HCHA)
pI_HC_HCHA_data21 <- data.frame(ID= pI_HC_HCHA_name1$.,seq=pI_HC_HCHA_data1$x, ID2= substr( pI_HC_HCHA_name1$., 1, 14))
pI_HC_HCHA_seq1=data.frame(pI_HC_HCHA_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HC_HCHA_data21$seq),1)) 
p_down_HC_HCHA_heatmap1<- p_down_HC_HCHA_heatmap[,2:21]
p_down_HC_HCHA_heatmap1$pI<- pI_HC_HCHA_seq1$PI1
p_down_HC_HCHA_heatmap1$ann<- p_down_HC_HCHA_heatmap$Nr.annotation
p_down_HC_HCHA_heatmap1$Protein<- p_down_HC_HCHA_heatmap$Protein
p_down_HC_HCHA_heatmap2 <- p_down_HC_HCHA_heatmap1[order(p_down_HC_HCHA_heatmap1$pI), ]
p_down_HC_HCHA_heatmap2_mRNA <-  normalized_hg[p_down_HC_HCHA_heatmap2$Protein,] 
#protein heatmap
down_HC_HCHA_heat_p=pheatmap(
  p_down_HC_HCHA_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
down_HC_HCHA_heat_mRNA=pheatmap(
  p_down_HC_HCHA_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
down_HC_HCHA_heat_p_bar <- ggplot(p_down_HC_HCHA_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_minimal() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线


#HCHA down 
HCHA_sequences <- entrez_fetch(db = "protein", id = p_down_HCHA_heatmap$Nr.ID, rettype = "fasta")
write(HCHA_sequences, file="HCHA_sequence.fasta", sep="\n")
pI_HCHA <- readBStringSet("HCHA_sequence.fasta", format = "fasta", use.names = TRUE)
pI_HCHA_name1 <- names(pI_HCHA) %>% as.data.frame()
pI_HCHA_data1 <- as.data.frame(pI_HCHA)
pI_HCHA_data21 <- data.frame(ID= pI_HCHA_name1$.,seq=pI_HCHA_data1$x, ID2= substr( pI_HCHA_name1$., 1, 14))
pI_HCHA_seq1=data.frame(pI_HCHA_data21) %>%
  mutate("PI1" = round(pI(seq = pI_HCHA_data21$seq),1)) 
p_down_HCHA_heatmap1<- p_down_HCHA_heatmap[,2:21]
p_down_HCHA_heatmap1$pI<- pI_HCHA_seq1$PI1
p_down_HCHA_heatmap1$ann<- p_down_HCHA_heatmap$Nr.annotation
p_down_HCHA_heatmap1$Protein<- p_down_HCHA_heatmap$Protein
p_down_HCHA_heatmap2 <- p_down_HCHA_heatmap1[order(p_down_HCHA_heatmap1$pI), ]
p_down_HCHA_heatmap2_mRNA <-  normalized_hg[p_down_HCHA_heatmap2$Protein,] 
#protein heatmap
down_HCHA_heat_p=pheatmap(
  p_down_HCHA_heatmap2[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#mRNA heatmap
down_HCHA_heat_mRNA=pheatmap(
  p_down_HCHA_heatmap2_mRNA[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = FALSE,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
down_HCHA_heat_p_bar <- ggplot(p_down_HCHA_heatmap2, aes(x = pI, y =  reorder(ann, -pI) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue", position = "dodge") +
  theme_classic() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")   # 在x=7处画一条虚线
all_plots <- plot_grid(down_HC_heat_mRNA$gtable, down_HC_heat_p$gtable, 
                       down_HC_HCHA_heat_mRNA$gtable, down_HC_HCHA_heat_p$gtable, 
                       down_HCHA_heat_mRNA$gtable, down_HCHA_heat_p$gtable, ncol = 2, align = "hv"  # "hv" stands for horizontal and vertical alignment
)
bar_plot<-plot_grid(down_HC_heat_p_bar,down_HC_HCHA_heat_p_bar,down_HCHA_heat_p_bar,ncol = 1,align = "hv")
ggsave("pI_down_combined_heatmaps.pdf", all_plots, width =6, height = 15)
ggsave("pI_down_combined_bar.pdf", bar_plot, width = 9, height = 15)


save.image(file = "D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
