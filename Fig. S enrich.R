# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c('ggdendro'))


library(AnnotationHub)
library(biomaRt)
library(clusterProfiler)
library("topGO")
library("Rgraphviz")
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
library(WGCNA)
library(ggcorrplot)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/enrich_protein')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/gene and protein correlation/corhghp1.RData")
traits <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/traits.csv', header = T, row.names=1) 
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")
#protein GO------------
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/DEP WGCNA/protein_WGCNA.R.RData")

protein_data2=data.frame(All_Gene_KME[which(All_Gene_KME$col=='black'),1])
protein_data1=data.frame(All_Gene_KME[which(All_Gene_KME$col=='greenyellow'),1])
#bprotein GO-of balck module----
proteinlist2=protein_data2[,1]
ego_black_BP <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "BP", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_black_GO <-data.frame(ego_black_BP@result) 

pdf("black_GO_protein_barplot.pdf",width=4,height=3) 
ggplot(data_all_black_GO [1:10,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("black_GO_protein_net.pdf",width=4,height=3) 
plotGOgraph(ego_black_BP,useInfo ="def")#GO-BP功能网络图
dev.off() 
#protein  kegg of balck module------------
#mRNA WGCNA KEGG PLOT-------
load("kegg_info.RData")
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")
gene_data1=data.frame(All_Gene_KME[which(All_Gene_KME$col=='black'),1])
gene_data2=data.frame(All_Gene_KME[which(All_Gene_KME$col=='brown'),1])


pathway2gene <- AnnotationDbi::select(Porites.orgdb,
                                      keys = keys(Porites.orgdb),
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

proteinlist2=protein_data2[,1]
ekp_black2 <- enricher(proteinlist2, 
                       TERM2GENE = pathway2gene, 
                       TERM2NAME = pathway2name, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff =1,
                       pAdjustMethod = "BH",
                       minGSSize =5)
ekp_black_1 <- ekp_black2@result[ekp_black2@result$pvalue <0.05&ekp_black2@result$Count>4,]
# ekp_black_2 <- ekp_black_1[-grep(pattern="disease|Renal|cancer|cardio|Leuk|leuk",ekp_black_1$Description),]
#delete pathway about disease
data_black_kegg <- data.frame(ekp_black_1)
pdf("black_KEGG_protein_barplot.pdf",width=5,height=8) 
ggplot(ekp_black_1[1:15,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 

# protein GO of greenyellow module------------
proteinlist2 = protein_data1[,1]
ego_greenyellow_BP <- enrichGO(gene = proteinlist2,
                               keyType = "GID", 
                               OrgDb =Porites.orgdb,
                               ont = "BP", #ALL BP或MF或CC
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05) 
data_greenyellow_GO <-data.frame(ego_greenyellow_BP@result) 

pdf("greenyellow protein GO_barplot.pdf",width=4,height=3) 
ggplot(data_greenyellow_GO [1:10,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("greenyellow_GO_protein_net.pdf",width=4,height=3) 
plotGOgraph(ego_greenyellow_BP,useInfo ="def")#GO-BP功能网络图
dev.off() 



#protein kegg of greenyellow module-------------
proteinlist2 = protein_data1[,1]
ekp_greenyellow2 <- enricher(proteinlist2, 
                             TERM2GENE = pathway2gene, 
                             TERM2NAME = pathway2name, 
                             pvalueCutoff =0.05, 
                             qvalueCutoff =1,
                             pAdjustMethod = "BH",
                             minGSSize =5)
ekp_greenyellow_2 <- ekp_greenyellow2@result[ekp_greenyellow2@result$pvalue <0.05&ekp_greenyellow2@result$Count>4,]

#delete pathway about disease
pdf("greenyellow_protein_kegg_barplot.pdf",width=5,height=6) 
ggplot(ekp_greenyellow_2[1:15,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#fig4 carbohydrate metabolic process-----
carbohydrate<- data_greenyellow_GO %>% filter (Description =="carbohydrate metabolic process")

carbohydrate_rownames<-data.frame(strsplit(carbohydrate$geneID,'/'))
p_carbohydrate_heatmap<-protein_data_1[rownames(protein_data_1) %in% carbohydrate_rownames[,1], ]
For_cor_carbohydrate<- data.frame(t(p_carbohydrate_heatmap),CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)
cor_carbohydrate<- cor(For_cor_carbohydrate)
##默认绘图square

cor_carbohydrate <- subset(cor_carbohydrate,select = -c(CO2,HCO3,CO3,pH,alk))
p_carbohydrate_heatmap<-cor_carbohydrate[rownames(cor_carbohydrate) %in% c("CO2","CO2","CO3","pH", "alk" ), ]

cor_carbohydrate1<- data.frame(round(t(p_carbohydrate_heatmap),2))

cor_carbohydrate1<- subset(cor_carbohydrate1,abs(CO2)>0.5|abs(alk)>0.5)
cor_carbohydrate2<- cor_carbohydrate1

cor_carbohydrate1[abs(cor_carbohydrate1)<0.5]="-"
carbohydrate_rownames_all=datahgraw1[rownames(datahgraw1) %in% row.names(cor_carbohydrate2), ]
cor_carbohydrate_Nr<- carbohydrate_rownames_all$Nr.annotation


pdf("carbohydrate.heatmap.pdf", width =4, height =6)
labeledHeatmap(Matrix = cor_carbohydrate2,
               xLabels = colnames(cor_carbohydrate2),
               yLabels = carbohydrate_rownames_all$Nr.annotation,
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix =cor_carbohydrate1,
               setStdMargins = FALSE,
               cex.text = 0.5)
dev.off()

##fig4-----translation-----
# translation<- data_all_black_GO  %>% filter (Description =="translation")
translation_rownames=data.frame(rownames(datahgraw1[grep(pattern="translation",datahgraw1$GO.BiologicalProcess),]))
# translation_rownames<-data.frame(strsplit(translation$geneID,'/'))
protein_data_1<- as.data.frame(protein_data_1)
p_translation_heatmap<-protein_data_1[rownames(protein_data_1) %in% translation_rownames[,1], ]
For_cor_translation<- data.frame(t(p_translation_heatmap),CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)
cor_translation<- cor(For_cor_translation)
##默认绘图square
cor_translation <- subset(cor_translation,select = -c(CO2,HCO3,CO3,pH,alk))
p_translation_heatmap<-cor_translation[rownames(cor_translation) %in% c("CO2","HCO3","CO3","pH", "alk" ), ]

cor_translation1<- data.frame(round(t(p_translation_heatmap),2))
cor_translation1<- subset(cor_translation1,abs(CO2)>0.5|abs(alk)>0.5)
cor_translation1<- cor_translation1[-5,]
cor_translation2<- cor_translation1

cor_translation1[abs(cor_translation1)<0.5]="-"
translation_rownames_all=datahgraw1[rownames(datahgraw1) %in% row.names(cor_translation2), ]
cor_translation_Nr<- translation_rownames_all$Nr.annotation

pdf("translation.heatmap.pdf", width =4, height =3)
labeledHeatmap(Matrix = cor_translation2,
               xLabels = colnames(cor_translation2),
               yLabels = translation_rownames_all$Nr.annotation,
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix =cor_translation1,
               setStdMargins = FALSE,
               cex.text = 0.5)
dev.off()

# hexose metabolic process----

hexose=datahgraw1[grep(pattern="hexose metabolic process",datahgraw1$GO.BiologicalProcess),]

hexose_rownames=rownames(hexose)
protein_data_1<- as.data.frame(protein_data_1)
p_hexose_heatmap<-protein_data_1[rownames(protein_data_1) %in% hexose_rownames, ]
hexose_P<-hexose[rownames(hexose) %in% rownames(p_hexose_heatmap), ]
pdf("protein_hexose.heatmap.pdf", width =5*(1+max(nchar(hexose_P$Nr.annotation))/90), height =10*length(hexose_P$Nr.annotation)/38)
pheatmap(p_hexose_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE, labels_row=c(hexose_P$Nr.annotation))
dev.off()