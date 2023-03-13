
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

setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')

##DEG------
datahgraw1 <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/All_Unigene.advance.annotation.csv', header = T, row.names=1) 
datahgraw <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/all.sample.annot.xls.csv', header = T, row.names=1) 
datahg1<- merge(datahgraw[,1:20],datahgraw1,by="row.names")
rownames(datahg1)<- datahg1[,1]
datahg1<- datahg1[,-1]
datahg<- datahg1[,1:20]
datahg[,5]<- rowMeans(datahg1[,6:8])#as the 5th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg[,18]<- rowMeans(datahg1[,19:20])#as the 18th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg <-round(as.matrix(datahg))
head(datahg)
# annotation exctraction
h_annot <- datahgraw1[,1:37]

#####DEG  
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

###gene normalization data
rldhg <- rlogTransformation(ddshg)  ## 得到经过DESeq2软件normlization 
normalized_hg=assay(rldhg)
boxplot(normalized_hg, col = rainbow(ncol(normalized_hg)*1.2),main="expression value",las=2)

## negative and positive DEGs-----------
hgConCaOH <- results(ddshg, contrast = c("conditionshg","Ca(OH)2","Control"))
hgConCO2 <- results(ddshg, contrast = c("conditionshg","CO2","Control"))
hgConCO2CaOH <- results(ddshg, contrast = c("conditionshg","CO2+Ca(OH)2","Control"))
hgConCO2NaOH <- results(ddshg, contrast = c("conditionshg","CO2+NaOH","Control"))

# DEP-------------------
protein_data_r <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/protein_LFQ.csv', header = T) 

datahgraw <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/All_Unigene.advance.annotation.csv', header = T) 
datahgraw_n<- data.frame(datahgraw[,1])
colnames(datahgraw_n) <- c("ID")
h_annot <- datahgraw[,1:37]
##merge with host gene to find host protein

protein_data<-merge(protein_data_r,datahgraw_n, by.x='Protein',by.y='ID')
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
protein_data %>% group_by(Protein) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
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
# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Control")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.2))
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
# Column names of the results table
colnames(data_results)

##positive and negative DEPs----

CaOH2_up <- data_results[data_results$CaOH2_vs_Control_p.adj<0.05&data_results$CaOH2_vs_Control_ratio>=1.2,][grep("^CaOH2_vs|name", colnames(data_results))]
CaOH2_down <- data_results[data_results$CaOH2_vs_Control_p.adj<0.05&data_results$CaOH2_vs_Control_ratio<=-1.2,][grep("^CaOH2_vs|name", colnames(data_results))]

CO2_up <- data_results[data_results$CO2_vs_Control_p.adj<0.05&data_results$CO2_vs_Control_ratio>=1.2,][grep("^CO2_vs|name", colnames(data_results))]
CO2_down <- data_results[data_results$CO2_vs_Control_p.adj<0.05&data_results$CO2_vs_Control_ratio<=-1.2,][grep("^CO2_vs|name", colnames(data_results))]

CO2CaOH2_up <- data_results[data_results$CO2CaOH2_vs_Control_p.adj<0.05&data_results$CO2CaOH2_vs_Control_ratio>=1.2,][grep("^CO2CaOH2_vs|name", colnames(data_results))]
CO2CaOH2_down <- data_results[data_results$CO2CaOH2_vs_Control_p.adj<0.05&data_results$CO2CaOH2_vs_Control_ratio<=-1.2,][grep("^CO2CaOH2_vs|name", colnames(data_results))]

CO2NaOH_up <- data_results[data_results$CO2NaOH_vs_Control_p.adj<0.05&data_results$CO2NaOH_vs_Control_ratio>=1.2,][grep("^CO2NaOH_vs|name", colnames(data_results))]
CO2NaOH_down <- data_results[data_results$CO2NaOH_vs_Control_p.adj<0.05&data_results$CO2NaOH_vs_Control_ratio<=-1.2,][grep("^CO2NaOH_vs|name", colnames(data_results))]

CaOH2_up <- data_results[data_results$CaOH2_vs_Control_p.adj<0.05&data_results$CaOH2_vs_Control_ratio>=1.2,][grep("^CaOH2_vs|name", colnames(data_results))]
CaOH2_down <- data_results[data_results$CaOH2_vs_Control_p.adj<0.05&data_results$CaOH2_vs_Control_ratio<=-1.2,][grep("^CaOH2_vs|name", colnames(data_results))]

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
                      show_percentage =T,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))

plot_list[[2]]=ggvenn(CO2_VENN,
                      show_percentage =T,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))

plot_list[[3]]=ggvenn(CO2CaOH2_VENN,
                      show_percentage =T,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb"),
                      set_name_color = c("#ff0000","#4a9b83"))
plot_list[[4]]=ggvenn(CO2NaOH_VENN,
                      show_percentage =T,
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



#VENN plot:protein CO2 vs AE ------- 


#CO2
# positive
genelist_CO2_positive= genelist_CO2_p[genelist_CO2_p$pvalue<0.05&genelist_CO2_p$FC> 1.2,]
# negative
genelist_CO2_negative= genelist_CO2_p[genelist_CO2_p$pvalue<0.05&genelist_CO2_p$FC< (-1.2),]

#CO2CaOH2
# positive
genelist_CO2CaOH2_positive= genelist_CO2CaOH2_p[genelist_CO2CaOH2_p$pvalue<0.05&genelist_CO2CaOH2_p$FC> 1.2,]
# negative
genelist_CO2CaOH2_negative= genelist_CO2CaOH2_p[genelist_CO2CaOH2_p$pvalue<0.05&genelist_CO2CaOH2_p$FC< (-1.2),]


##positive
##CO2NaOH
# positive
genelist_CO2NaOH_positive= genelist_CO2NaOH_p[genelist_CO2NaOH_p$pvalue<0.05&genelist_CO2NaOH_p$FC> 1.2,]
# negative
genelist_CO2NaOH_negative= genelist_CO2NaOH_p[genelist_CO2NaOH_p$pvalue<0.05&genelist_CO2NaOH_p$FC< (-1.2),]

#reverse_up_gene 
genelist_CO2p<-genelist_CO2_positive$GID
genelist_CO2CaOHp<-genelist_CO2CaOH2_positive$GID
genelist_CO2NaOHp<-genelist_CO2NaOH_positive$GID
gene_reverse<-list(`CO2`=genelist_CO2p, 
                   `CO2CaOH`= genelist_CO2CaOHp,
                   `CO2NaOH`= genelist_CO2NaOHp)

#reverse_down_gene
genelist_CO2d<-genelist_CO2_negative$GID
genelist_CO2CaOHd<-genelist_CO2CaOH2_negative$GID
genelist_CO2NaOHd<-genelist_CO2NaOH_negative$GID
gene_reverse_d<-list(`CO2`=genelist_CO2d, 
                     `CO2CaOH`= genelist_CO2CaOHd,
                     `CO2NaOH`= genelist_CO2NaOHd)


inter_down <- get.venn.partitions(gene_reverse_d)

#venn.partitions and enrich-----
library(AnnotationHub)
library(biomaRt)
library(clusterProfiler)
library("topGO")
library("Rgraphviz")
library(KEGGREST)
library(AnnotationForge)
library(splitstackshape)
library(tidyverse)
library(jsonlite)
library(enrichplot)
library(ggstance)
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")
#up 
inter_up <- get.venn.partitions(gene_reverse)
inter_up_reversible <- unlist(inter_up[7,5])
ego_inter_up_reversible<- enrichGO(gene = inter_up_reversible,
                                   keyType = "GID", 
                                   OrgDb =Porites.orgdb,
                                   ont = "BP", #ALL BP或MF或CC
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  =0.05) 

inter_up_by <- unlist(inter_up[2,5])
ego_inter_up_by <- enrichGO(gene = inter_up_by,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "BP", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  =0.05) 

inter_up_not<- unlist(inter_up[1,5])
ego_inter_up_not <- enrichGO(gene = inter_up_not,
                             keyType = "GID", 
                             OrgDb =Porites.orgdb,
                             ont = "BP", #ALL BP或MF或CC
                             pAdjustMethod = "BH",
                             pvalueCutoff  =0.05) 

#down
inter_down <- get.venn.partitions(gene_reverse_d)
inter_down_reversible <- unlist(inter_down[7,5])
ego_inter_down_reversible<- enrichGO(gene = inter_down_reversible,
                                     keyType = "GID", 
                                     OrgDb =Porites.orgdb,
                                     ont = "BP", #ALL BP或MF或CC
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  =0.05) 
inter_down_by <- unlist(inter_down[2,5])
ego_inter_down_by <- enrichGO(gene = inter_down_by,
                              keyType = "GID", 
                              OrgDb =Porites.orgdb,
                              ont = "BP", #ALL BP或MF或CC
                              pAdjustMethod = "BH",
                              pvalueCutoff  =0.05) 


inter_down_not<- unlist(inter_down[1,5])
ego_inter_down_not <- enrichGO(gene = inter_down_not,
                               keyType = "GID", 
                               OrgDb =Porites.orgdb,
                               ont = "BP", #ALL BP或MF或CC
                               pAdjustMethod = "BH",
                               pvalueCutoff  =0.05) 
##Fig 2B -up ven barplot------- 
plot_list = list() 

plot_list[[1]]=ggvenn(gene_reverse,
                      show_percentage =T,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
                      set_name_color = c("#ff0000","#4a9b83","#1d6295"))


plot_list[[3]]=ggplot(ego_inter_up_by@result[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + labs(title="side effects")+ 
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 

plot_list[[4]]=ggplot(ego_inter_up_not@result[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + labs(title="irreversible")+
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 



pdf("Fig 2B GO_up_barplot.pdf",width=40,height=10) 
pushViewport(viewport(layout = grid.layout(1,4)))
print(plot_list[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col =1))
print(plot_list[[2]], vp = viewport(layout.pos.row = 1, layout.pos.col =2))
print(plot_list[[3]], vp = viewport(layout.pos.row = 1, layout.pos.col =3))
print(plot_list[[4]], vp = viewport(layout.pos.row = 1, layout.pos.col =4))

dev.off() 


#Fig 2C-down ven barplot------- 

plot_list[[5]]=ggvenn(gene_reverse_d,
                      show_percentage =T,
                      stroke_color = "NA",
                      fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
                      set_name_color = c("#ff0000","#4a9b83","#1d6295"))
plot_list[[6]]=ggplot(ego_inter_down_reversible@result[1:5,], aes(-log(pvalue),fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + labs(title="reversible")+
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 

plot_list[[7]]=ggplot(ego_inter_down_by@result[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity')+ labs(title="side effects")+ 
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 

plot_list[[8]]=ggplot(ego_inter_down_not@result[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity')+ labs(title="irreversible")+
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 

pdf("Fig 2C GO_down_barplot.pdf",width=40,height=10) 
pushViewport(viewport(layout = grid.layout(1,4)))
print(plot_list[[5]], vp = viewport(layout.pos.row = 1, layout.pos.col =1))
print(plot_list[[6]], vp = viewport(layout.pos.row = 1, layout.pos.col =2))
print(plot_list[[7]], vp = viewport(layout.pos.row = 1, layout.pos.col =3))
print(plot_list[[8]], vp = viewport(layout.pos.row = 1, layout.pos.col =4))
dev.off() 
save.image("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")