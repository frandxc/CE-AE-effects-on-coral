setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/gene and protein correlation')
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
require(cowplot)
library(gridExtra)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(amap)
library(ggpubr)
library(purrr)
library(pheatmap) 
library(limma)
library(sva)
library(bladderbatch)
library(stringr)
library(psych)
load('corhghp1.RData')


# mRNA phosphoglycerate kinase-----
phosphoglycerate=datahgraw1[grep(pattern="phosphoglycerate kinase",datahgraw1$Nr.annotation),]
phosphoglycerate_rownames=rownames(phosphoglycerate)
normalized_hg_d<- as.data.frame(normalized_hg)
phosphoglycerate_heatmap<-normalized_hg_d[rownames(normalized_hg_d) %in% phosphoglycerate_rownames, ]
phosphoglycerate_heatmap1<-phosphoglycerate_heatmap[rowSums(phosphoglycerate_heatmap[,1:20])>150,]


phosphoglycerate_heatmap_g<-phosphoglycerate_heatmap[rownames(phosphoglycerate_heatmap) %in% rownames(phosphoglycerate_heatmap1), ]

phosphoglycerate_MEAN1=colMeans(phosphoglycerate_heatmap1)
phosphoglycerate_data=data.frame(conditionshg,phosphoglycerate_MEAN1)
phosphoglycerate_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pheatmap(phosphoglycerate_heatmap_g,scale='row',border_color=NA,cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE, labels_row=c(phosphoglycerate$Nr.annotation))

# PROTEIN  phosphoglycerate-----
p_phosphoglycerate_heatmap<-protein_data_1[rownames(protein_data_1) %in% phosphoglycerate_rownames, ]

p_phosphoglycerate_heatmap1<-p_phosphoglycerate_heatmap[rowSums(p_phosphoglycerate_heatmap[,1:20])>1,]

p_phosphoglycerate_heatmap_g<-p_phosphoglycerate_heatmap[rownames(p_phosphoglycerate_heatmap) %in% rownames(p_phosphoglycerate_heatmap1), ]

p_phosphoglycerate_MEAN1=colMeans(p_phosphoglycerate_heatmap1)

p_phosphoglycerate_data=data.frame(conditionshg,p_phosphoglycerate_MEAN1)
colnames(p_phosphoglycerate_data)
p_phosphoglycerate_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pheatmap(p_phosphoglycerate_heatmap_g,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE, labels_row=c(phosphoglycerate$Nr.annotation))




#tRNA ligase---- 
tRNA_ligase=datahgraw1[grep(pattern="tRNA ligase",datahgraw1$Nr.annotation),]
tRNA_ligase_rownames=rownames(tRNA_ligase)
normalized_hg_d<- as.data.frame(normalized_hg)
tRNA_ligase_heatmap<-normalized_hg_d[rownames(normalized_hg_d) %in% tRNA_ligase_rownames, ]
tRNA_ligase_heatmap1<-tRNA_ligase_heatmap[rowSums(tRNA_ligase_heatmap[,1:20])>1,]


tRNA_ligase_heatmap_g<-tRNA_ligase_heatmap[rownames(tRNA_ligase_heatmap) %in% rownames(tRNA_ligase_heatmap1), ]

tRNA_ligase_MEAN1=colMeans(tRNA_ligase_heatmap1)

tRNA_ligase_data=data.frame(conditionshg,tRNA_ligase_MEAN1)
tRNA_ligase_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pheatmap(tRNA_ligase_heatmap_g,scale='row',border_color=NA,cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE, labels_row=c(tRNA_ligase$Nr.annotation))

# protein tRNA_ligase -----
p_tRNA_ligase_heatmap<-protein_data_1[rownames(protein_data_1) %in% tRNA_ligase_rownames, ]
p_tRNA_ligase_heatmap1<-p_tRNA_ligase_heatmap
p_tRNA_ligase_heatmap_g<-tRNA_ligase[rownames(p_tRNA_ligase_heatmap), ]

p_tRNA_ligase_MEAN1=colMeans(p_tRNA_ligase_heatmap1)


p_tRNA_ligase_data=data.frame(conditionshg,p_tRNA_ligase_MEAN1)
colnames(p_tRNA_ligase_data)
p_tRNA_ligase_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pheatmap(p_tRNA_ligase_heatmap1,scale='row',border_color=NA,cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE, labels_row=c(p_tRNA_ligase_heatmap_g$Nr.annotation))
#PCA------
traits <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/traits.csv', header = T, row.names=1) 
CAS_data=data.frame(conditionshg,traits$CAS)
CAS_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
#phosphoglycerate---- 
phosphoglycerate_data=data.frame(conditionshg,traits$phosphoglycerate)
phosphoglycerate_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#tRNAligase---- 
tRNAligase_data=data.frame(conditionshg,traits$tRNAligase)
tRNAligase_data$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

library(factoextra)
citation("factoextra")

#phosphoglyceratepca
traitsphosphoglycerate=data.frame(traits$phosphoglycerate,phosphoglycerate_MEAN1,p_phosphoglycerate_MEAN1)
traitsphosphoglycerate.pca <- prcomp(traitsphosphoglycerate, scale. = TRUE)


#tRNAligase pca
traitstRNAligase=data.frame(traits$tRNAligase,tRNA_ligase_MEAN1,p_tRNA_ligase_MEAN1)
traitstRNAligase.pca <- prcomp(traitstRNAligase, scale. = TRUE)


#fig for supplementary------
compaired <- list(c("Control", "Ca(OH)2"), 
                  c("Control","CO2"), 
                  c("Control","CO2+Ca(OH)2"),
                  c("Control","CO2+NaOH"))
fig3=ggplot(data = phosphoglycerate_data,aes(x=conditionshg ,y=phosphoglycerate_MEAN1),ordered=TRUE)+labs(title = "phosphoglycerate")+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)

fig4=ggplot(data = p_phosphoglycerate_data,aes(x=conditionshg ,y=p_phosphoglycerate_MEAN1),ordered=TRUE)+labs(title = "phosphoglycerate")+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)


fig5=ggplot(data = tRNA_ligase_data,aes(x=conditionshg ,y=tRNA_ligase_MEAN1),ordered=TRUE)+labs(title = "tRNA_ligase")+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)

fig6=ggplot(data = p_tRNA_ligase_data,aes(x=conditionshg ,y=p_tRNA_ligase_MEAN1),ordered=TRUE)+labs(title = "tRNA_ligase")+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)



# fig for text-----

fig8=ggplot(data = phosphoglycerate_data,aes(x=conditionshg ,y=traits.phosphoglycerate),ordered=TRUE)+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)

fig9=ggplot(data = tRNAligase_data,aes(x=conditionshg ,y=traits.tRNAligase),ordered=TRUE)+
  geom_boxplot()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)
#pca

fig11=fviz_pca_var(traitsphosphoglycerate.pca ,repel = TRUE, col.circle ="transparent",  col.var = "cos2",gradient.cols = c("black", "blue", "red"))+  theme_classic()+ theme(legend.position = 'none')

fig12=fviz_pca_var(traitstRNAligase.pca,repel = TRUE, col.circle ="transparent",  col.var = "cos2",gradient.cols = c("black", "blue", "red"))+  theme_classic()+ theme(legend.position = 'none')


library(FactoMineR)
require(cowplot)
library(gridExtra)
library(grid)
pdf("enzyme_plots.pdf",  width =10, height =8) #paper= "a4"
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,4)))
print(fig3, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(fig4, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(fig8, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(fig11, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(fig5, vp = viewport(layout.pos.row =2, layout.pos.col = 1))
print(fig6, vp = viewport(layout.pos.row =2, layout.pos.col = 2))
print(fig9, vp = viewport(layout.pos.row =2, layout.pos.col = 3))
print(fig12, vp = viewport(layout.pos.row =2, layout.pos.col = 4))

dev.off()  



