# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c('ggdendro'))

# load package------------
# load library------------
library(AnnotationHub)
library(biomaRt)
library(clusterProfiler)
library("topGO")
library("Rgraphviz")
library("pathview")
library('stringr')
library(dplyr)
library(KEGGREST)
library(AnnotationForge)
library(splitstackshape)
library(tidyverse)
library(jsonlite)
library(purrr)
library(RCurl)
library(stringr)
library(enrichplot)
library(ggnewscale)
library(forcats)
library(ggstance)
library(GOplot)
library(ggrepel)
library(ggvenn)
library("futile.logger")
library("VennDiagram")
library('ggplot2')
library(ggpubr)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")

load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")

mediancenter <- function(x, refrows=NULL) {
  if (!is.null(refrows)) {
    gm <- apply(x[refrows,], 2, median, na.rm = T)
    x <- sweep(x, 2, gm, `-`)+median(x[refrows,], na.rm=T)
  } else {
    gm <- apply(x, 2, median, na.rm=T)
    x <- sweep(x, 2, gm, `-`)+median(x, na.rm=T)
  }
  return(x)
}
# gene TPM normalization------
datahg2<- datahg1
datahg2[,5]<- rowMeans(datahg2[,6:8])#as the 5th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg2[,18]<- rowMeans(datahg2[,19:20])#as the 18th sample was an outlier, and it is replaced by the mean of 19 and 20th samples
datahg3 <-round(as.matrix(datahg2[,1:20]))
kb <-datahg2$geneLength/1000
datahg_n1<-datahg3[,1:20]/kb
tpm_g <- t(t(datahg_n1)/colSums(datahg_n1) * 1000000)
# log2 transformed
tpm_g1 <- data.frame(tpm_g,rowMin(tpm_g))
names(tpm_g1)
tpm_g2 <-tpm_g1 %>%filter(rowMin.tpm_g.>1) # to avoid the inf value
tpm_g_log <-log2(tpm_g+1)
tpm_g_log_cen <- data.frame(mediancenter(tpm_g_log[,1:20]))

#protein mediancenter transformed -------------------
#####protein data filter------

protein_data_1 <- data_imp@assays@data@listData[[1]]
conditionshg <- c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2Ca(OH)2",4),rep("CO2NaOH",4))
replicates=4
datahp_log_cen <- mediancenter(protein_data_1)

# merge gene and protein----------------
hghp_merge<- merge(tpm_g_log_cen,datahp_log_cen,by="row.names")
row.names(hghp_merge) <- hghp_merge$Row.names
hghp_merge_f <- hghp_merge[,-1] 
hghp_merge_f <- hghp_merge[,-1] 
hghp_merge_f[hghp_merge_f=="-Inf"] <-0 

hghp_annot_merge<- merge(hghp_merge,h_annot,by=1)
row.names(hghp_annot_merge) <- hghp_annot_merge$Row.names
hghp_annot_merge_f1 <- hghp_annot_merge[,-1] 


# mRNA/protein correlation between treatments------------
x <- hghp_merge_f[,1:40]
corhghp_each <- matrix(data=0, ncol=1, nrow = length(hghp_merge_f[1,])/2)
# hgnor_t<-  t(hghp_merge_f[,1:20])
# hpnor_t<-  t(hghp_merge_f[,21:40])

for (j in 1:20) 
{
  corhghp_each[j,]<-  cor(hghp_merge_f[,j],hghp_merge_f[,j+20])
}
corhghp_each_m <- matrix(corhghp_each, ncol = 4, byrow = TRUE)
averge_gp<-apply(corhghp_each_m, 1, mean)
sd_gp<-apply(corhghp_each_m,1, sd)  
ave_sd<-cbind(averge_gp,sd_gp)

colnames(ave_sd)<-c('average','sd')
rownames(ave_sd)<-c('Control','Ca(OH)2','CO2','CO2+Ca(OH)2','CO2+NaOH')

#Fig 3 A SCATTER mRNA_protein-----
hghp_merge_f_c<-  c(hghp_merge_f[,1])
hghp_merge_f_c_p<-  c(hghp_merge_f[,1+20])
cor(hghp_merge_f_c,hghp_merge_f_c_p)
scatterdata_gp<-  data.frame(hghp_merge_f_c,hghp_merge_f_c_p)
colnames(scatterdata_gp)<- c('gene','protein')
scatterdata_gp <-  scatterdata_gp%>%filter(gene>2)  
pdf(file="Fig3 A_CATTER mRNA_protein.pdf")
ggplot(scatterdata_gp,aes(gene,protein),)+geom_point()+  theme_classic()
dev.off()


#Fig 3Btotal mRNA----
compaired <- list(c("Con", "CaOH"), 
                  c("Con","CO2"), 
                  c("Con","CO2CaOH"),
                  c("Con","CO2NaOH"))
conditionshg <- factor(c(rep("Con",4),rep( "CaOH",4),rep("CO2",4),rep("CO2CaOH",4),rep("CO2NaOH",4)))

mean_mRNA_Protein=data.frame(colMeans(hghp_merge_f))
mean_data_mRNA <-data.frame(conditionshg,mRNA=mean_mRNA_Protein[1:20,] )
mean_data_mRNA$conditionshg <- factor(mean_data_mRNA$conditionshg, levels=c("Con", "CaOH","CO2","CO2CaOH","CO2NaOH"), ordered=TRUE)

pdf(file="Fig.3B_mean_mRNA_total.pdf")
mean_data_mRNA%>%  ggplot(aes(conditionshg,mRNA,color=conditionshg),palette = "jco", ordered=TRUE)+
  geom_boxplot()+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T,test = t.test) 
dev.off()

#Fig 3Ctotal protein----
mean_mRNA_Protein=data.frame(colMeans(hghp_merge_f))
mean_data_protein<-data.frame(conditionshg,protein=mean_mRNA_Protein[21:40,] )
mean_data_protein$conditionshg <- factor(mean_data_protein$conditionshg, levels=c("Con", "CaOH","CO2","CO2CaOH","CO2NaOH"), ordered=TRUE)

pdf(file="Fig.3C mean_protein_total.pdf")
mean_data_protein%>%  ggplot(aes(conditionshg,protein,color=conditionshg),palette = "jco", ordered=TRUE)+
  geom_boxplot()+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T,test = t.test) 

dev.off()
#Fig 3Dcorrelation_mRNA_protein_traits----
traits <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/traits.csv', header = T, row.names=1) 
For_cor_protein_mRNA<- data.frame(mRNA=mean_data_mRNA$mRNA,protein=mean_data_protein$protein,CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)
cor_protein_mRNA<- cor(For_cor_protein_mRNA, use = "p")
library(ggcorrplot)
library(WGCNA)
cor_protein_mRNA1<- data.frame(round(cor_protein_mRNA[1:2,3:7],3))
#p value of correlation
nSamples = nrow(For_cor_protein_mRNA);
cor_Pvalue = corPvalueStudent(cor_protein_mRNA, nSamples)
#cor_protein_mRNA1<- subset(cor_protein_mRNA1,abs(CO2)>0.5|abs(alk)>0.5)
cor_protein_mRNA2<- cor_protein_mRNA1
cor_protein_mRNA1[abs(cor_protein_mRNA1)<0.5]="-"
pdf("Fig.3D protein_mRNA.heatmap.pdf", width =4, height =4)
labeledHeatmap(Matrix = t(cor_protein_mRNA2),
               xLabels = rownames(cor_protein_mRNA2),
               yLabels = colnames(cor_protein_mRNA2),
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix = t(cor_protein_mRNA1),
               setStdMargins = FALSE,
               cex.text = 0.5)
dev.off()
#FIG 3E g&p correlation within genes------------
hgnor_t<-  t(hghp_merge_f[,1:20])
hpnor_t<-  t(hghp_merge_f[,21:40])
corhghp_gene <- matrix(data=0, nrow =length(hghp_merge_f[,1]), ncol = 1)
for (i in 1:length(hghp_merge_f[,1])) 
{
  corhghp_gene[i,]<-  cor(hgnor_t[,i],hpnor_t[,i])
}
rownames(corhghp_gene) <-hghp_merge$Row.names

#g&p plot cor betw genes---
pdf(file="Fig.3E HIST_HPHG_COR.pdf")
hist(corhghp_gene)
dev.off()

#enrich------
library(clusterProfiler)
library("topGO")
library("Rgraphviz")
library("pathview")
library('stringr')
library(dplyr)
library(KEGGREST)
library(AnnotationForge)
library(splitstackshape)
library(tidyverse)
library(jsonlite)
library(purrr)
library(RCurl)
library(stringr)
library(enrichplot)
library(ggnewscale)
library(forcats)
library(ggstance)
library(GOplot)
library(ggrepel)
library(ggvenn)
library("futile.logger")
library("VennDiagram")

Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/make orgdb/org.P_pukoensis.eg.sqlite")

#positive cor-----
cor_hgh_POS=subset(corhghp_gene,corhghp_gene[,1]>0.5) 
cor_hgh_NEG=subset(corhghp_gene,corhghp_gene[,1]<(-0.5)) 
cor_hgh_POS_list=rownames(data.frame(cor_hgh_POS))
cor_hgh_NEG_list=rownames(data.frame(cor_hgh_NEG))

cor_hgh_POS_list_ego <- enrichGO(gene =cor_hgh_POS_list,
                                 keyType = "GID", 
                                 OrgDb =Porites.orgdb,
                                 ont = "BP", #ALL BP或MF或CC
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05) 

cor_POS_GO <- cor_hgh_POS_list_ego@result[cor_hgh_POS_list_ego@result$pvalue <0.05&cor_hgh_POS_list_ego@result$Count>4,]

p1=ggplot(cor_POS_GO, aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9),text= element_text(size=9))+coord_fixed(ratio=4) 

#negative cor-----
cor_hgh_NEG_list_ego <- enrichGO(gene =cor_hgh_NEG_list,
                                 keyType = "GID", 
                                 OrgDb =Porites.orgdb,
                                 ont = "BP", #ALL BP或MF或CC
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05) 

cor_NEG_GO <- cor_hgh_NEG_list_ego@result[cor_hgh_NEG_list_ego@result$pvalue <0.05&cor_hgh_NEG_list_ego@result$Count>4,]

p2=ggplot(cor_NEG_GO, aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue), showCategory=10) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red') + 
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 

pdf("Fig.3FG_neg and posi_go.pdf",width=8,height=8) 
par(mfrow=c(1,2))
print(p1)
print(p2)
dev.off() 

#protein POSITIVE KEGG PLOT-------
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/kegg_info.RData")

Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")
pathway2gene <- AnnotationDbi::select(Porites.orgdb,
                                      keys = keys(Porites.orgdb),
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

ekp_pos2 <- enricher(cor_hgh_POS_list, 
                     TERM2GENE = pathway2gene, 
                     TERM2NAME = pathway2name, 
                     pvalueCutoff =0.05, 
                     qvalueCutoff =1,
                     pAdjustMethod = "BH",
                     minGSSize =5)
ekp_pos_1 <- ekp_pos2@result[ekp_pos2@result$pvalue <0.05&ekp_pos2@result$Count>4,]
# ekp_pos_2 <- ekp_pos_1[-grep(pattern="disease|Renal|cancer|cardio|Leuk|leuk",ekp_pos_1$Description),]
#delete pathway about disease
data_pos_kegg <- data.frame(ekp_pos_1)
pdf("pos_KEGG_protein_barplot.pdf",width=5,height=4) 
ggplot(ekp_pos_1[1:15,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9),text= element_text(size=9))+coord_fixed(ratio=4) 
dev.off() 

#protein NEGATIVE KEGG PLOT-------

ekp_NEG2 <- enricher(cor_hgh_NEG_list, 
                     TERM2GENE = pathway2gene, 
                     TERM2NAME = pathway2name, 
                     pvalueCutoff =0.05, 
                     qvalueCutoff =1,
                     pAdjustMethod = "BH",
                     minGSSize =5)
ekp_NEG_1 <- ekp_NEG2@result[ekp_NEG2@result$pvalue <0.05&ekp_pos2@result$Count>4,]
# ekp_NEG_2 <- ekp_NEG_1[-grep(pattern="disease|Renal|cancer|cardio|Leuk|leuk",ekp_NEG_1$Description),]
#delete pathway about disease
data_NEG_kegg <- data.frame(ekp_NEG_1)
pdf("NEG_KEGG_protein_barplot.pdf",width=5,height=8) 
ggplot(ekp_NEG_1[1:15,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9),text= element_text(size=9))+coord_fixed(ratio=4) 
dev.off() 

save.image("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig3.RData")
