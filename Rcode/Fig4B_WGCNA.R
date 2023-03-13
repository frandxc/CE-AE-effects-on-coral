#BiocManager::install("ggord")
#Protein and transcript datasets were median centred to the overall median of the respective dataset.
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

library("RColorBrewer")
library("amap")
library("ggplot2")
library("BiocParallel")
library(ggpubr)
library(purrr)
library(pheatmap) 
library(dplyr)
library(patchwork)
require(cowplot)
library(gridExtra)
library(gplots)
library(amap)
library(limma)
library(cluster)
library(ggfortify)
library(WGCNA)

setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig3.RData")
# annotation exctraction
# annotation exctraction--
#PTR:ratio of protein:transcripts-----
frorPTR<-  hghp_merge_f[-which(apply(hghp_merge_f,1,min)<1),]
proteinforratio<-  frorPTR[,21:40]
geneforratio<-  frorPTR[,1:20]
PTR<- proteinforratio/geneforratio

PTR_data_1 <- as.matrix(PTR)
conditionshg <- c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2Ca(OH)2",4),rep("CO2NaOH",4))
replicates=4
datahp_log_cen <- mediancenter(PTR_data_1)

cont_mean<-rowMeans(PTR_data_1[,1:4])
cont_mean_g<-rowMeans(geneforratio[,1:4])
cont_mean_p<-rowMeans(proteinforratio[,1:4])

co2_mean<-rowMeans(PTR_data_1[,9:12])
co2_mean_g<-rowMeans(geneforratio[,9:12])
co2_mean_p<-rowMeans(proteinforratio[,9:12])

ratio_mean<-data.frame(sort(co2_mean/cont_mean,decreasing=T))
list1<-rownames(ratio_mean)[1:100]
list2<-rownames(ratio_mean)[1528:1628]

plot(cont_mean_g,cont_mean)
plot(cont_mean_p,cont_mean)
plot(co2_mean,co2_mean_g)
plot(co2_mean,co2_mean_p)



#PTR WGCNA--------------------

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
#######input data
# hpraw <- assay(rldhg)
hpraw1 <- datahp_log_cen 
# hpraw1[hpraw1==	-Inf] <- 0
hpraw2 <- hpraw1[order(rowSums(hpraw1), decreasing=T), ]

datExpr0 = as.data.frame(t(hpraw2));
sampleTree = hclust(dist(datExpr0), method = "median");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# traitData------
traitData = read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/traits.csv', header = T, row.names=1);
dim(traitData)
names(traitData)
allTraits = traitData[,1:9];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data 
hpSamples = rownames(datExpr0);
traitRows = match(hpSamples, allTraits);
datTraits = allTraits;

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Sample dendrogram and trait heatmap.pdf", width = 10, height =6)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
#  Code chunk 9
dev.off()

#Step-by-step network construction and module detection
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, RsquaredCut = 0.85, powerVector = powers, verbose = 5)
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file='sft.pdf',w=10,h=7)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#  powerEstimate
softPower = sft$powerEstimate;
adjacency = adjacency(datExpr0, power = softPower);
#  Turn adjacency into topological overlap-
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#We like large modules, so we set the minimum module size relatively high:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 3, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#  Code chunk 7
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf(file='METree.pdf',w=10,h=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()
#the cut line into the dendrogram
MEDissThres = 0.8  #越小，颜色越多
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf(file = "DendroAndColors3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
#Relating modules to external clinical traits and identifying important genes

nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor1=signif(moduleTraitCor,1)
moduleTraitCor1[abs(moduleTraitCor1)<0.5]=' '

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue1=signif(moduleTraitPvalue,1)
moduleTraitPvalue1[moduleTraitPvalue1<=0.01]='*'
moduleTraitPvalue1[moduleTraitPvalue1>0.01]=' '

sizeGrWindow(20,20)

# genes in different module-------------------
blankdata <- data.frame()
write.table(blankdata,quote = F,row.names = T,file = "All_Gene_KME.txt",col.names = F)
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr0[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,sep='\t',row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}
All_Gene_KME<-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/All_Gene_KME.txt',sep='')

colnames(All_Gene_KME)<-c('id','val','col')
# wgcna.Module-trait.heatmap.pdf-------

#set y labs
MElength=data.frame(length(All_Gene_KME[,1]))
MElength[1,]=length(All_Gene_KME[which(All_Gene_KME$col=='blue'),1])
MElength[2,]=length(All_Gene_KME[which(All_Gene_KME$col=='black'),1])
MElength[3,]=length(All_Gene_KME[which(All_Gene_KME$col=='turquoise'),1])

MElength1=as.character(MElength[,1])
text1 =paste(substring(colnames(MEs),3),"\n",'n=',MElength1)
pdf("wgcna.PTR-trait.heatmap.pdf", width = 12, height =5)
textMatrix =  paste(moduleTraitCor1, "\n",
                    moduleTraitPvalue1, "", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 6, 1, 2));
labeledHeatmap(Matrix = moduleTraitCor[1:3,],
               xLabels = names(allTraits),
               yLabels = names(MEs[,1:3]),
               ySymbols = text1,
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 2,
               zlim = c(-1,1))
dev.off()

#FIG protein GO------------
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
library(ggcorrplot)
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/make orgdb/org.P_pukoensis.eg.sqlite")
#bprotein GO-of blue module----
protein_data2=data.frame(All_Gene_KME[which(All_Gene_KME$col=='blue'),1])

proteinlist2=protein_data2[,1]
ego_blue_BP <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "BP", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_blue_GO <-data.frame(ego_blue_BP@result) 
proteinlist2=protein_data2[,1]
ego_blue_CC <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "CC", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_blue_CC <-data.frame(ego_blue_CC@result) 

pdf("blue_GO_PTR_BP.pdf",width=8,height=3) 
ggplot(data_all_blue_GO [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="blue"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("blue_GO_PTR_CC.pdf",width=8,height=3) 
ggplot(data_all_blue_CC [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="blue"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#protein kegg of blue module-------------
ekp_blue2 <- enricher(proteinlist2, 
                             TERM2GENE = pathway2gene, 
                             TERM2NAME = pathway2name, 
                             pvalueCutoff =0.05, 
                             qvalueCutoff =1,
                             pAdjustMethod = "BH",
                             minGSSize =5)
ekp_blue_2 <- ekp_blue2@result[ekp_blue2@result$pvalue <0.05&ekp_blue2@result$Count>4,]

#delete pathway about disease
pdf("blue_protein_kegg_barplot.pdf",width=5,height=6) 
ggplot(ekp_blue_2[1:15,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#bprotein GO-of turquoise module----
protein_data2=data.frame(All_Gene_KME[which(All_Gene_KME$col=='turquoise'),1])
proteinlist2=protein_data2[,1]
ego_turquoise_BP <- enrichGO(gene = proteinlist2,
                          keyType = "GID", 
                          OrgDb =Porites.orgdb,
                          ont = "BP", #ALL BP或MF或CC
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05) 
data_all_turquoise_GO <-data.frame(ego_turquoise_BP@result) 
proteinlist2=protein_data2[,1]
ego_turquoise_CC <- enrichGO(gene = proteinlist2,
                          keyType = "GID", 
                          OrgDb =Porites.orgdb,
                          ont = "CC", #ALL BP或MF或CC
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05) 
data_all_turquoise_CC <-data.frame(ego_turquoise_CC@result) 

pdf("turquoise_GO_PTR_BP.pdf",width=8,height=3) 
ggplot(data_all_turquoise_GO [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="turquoise"),text= element_text(size=9,color="turquoise"))+coord_fixed(ratio=4) 
dev.off() 
pdf("turquoise_GO_PTR_CC.pdf",width=8,height=3) 
ggplot(data_all_turquoise_CC [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="turquoise"),text= element_text(size=9,color="turquoise"))+coord_fixed(ratio=4) 
dev.off() 