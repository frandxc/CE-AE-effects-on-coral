#BiocManager::install(c("str_split_fixed"),force = TRUE,lib = "D:/myR/library")
#Protein and transcript datasets were median centred to the overall median of the respective dataset.


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
# annotation exctraction
# annotation exctraction--
# datahgraw <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/gene_host.csv', header = T, row.names=1) 
# h_annot <- datahgraw[,21:30]
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
protein_data_1 <- data_imp@assays@data@listData[[1]]
protein_data_2 <-data.frame(protein_data_1[rownames(protein_data_1) %in% rownames(filtered_coral),])

conditionshg <- c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2Ca(OH)2",4),rep("CO2NaOH",4))
replicates=4
# mediancenter <- function(x, refrows=NULL) {
#   if (!is.null(refrows)) {
#     gm <- apply(x[refrows,], 2, median, na.rm = T)
#     x <- sweep(x, 2, gm, `-`)+median(x[refrows,], na.rm=T)
#   } else {
#     gm <- apply(x, 2, median, na.rm=T)
#     x <- sweep(x, 2, gm, `-`)+median(x, na.rm=T)
#   }
#   return(x)
# }

# 进行中心化处理
datahp_log_cen1 <- protein_data_1

datahp_log_cen <-data.frame(datahp_log_cen1[rownames(datahp_log_cen1) %in% rownames(filtered_coral),])
#protein WGCNA--------------------

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

# traitData------
traitData = read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/traits.csv', header = T, row.names=1);
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

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
pdf("Dendro.heatmap.pdf", width =5, height =4)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#  Code chunk 9

#Step-by-step network construction and module detection
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, RsquaredCut = 0.85, powerVector = powers, verbose = 5)
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#  powerEstimate
softPower = 7;
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
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit =1, pamRespectsDendro = FALSE,
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
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#the cut line into the dendrogram
MEDissThres = 0.6  #越小，颜色越多
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
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
MEs<- MEs[, 1:(ncol(MEs)-1)]
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor1=signif(moduleTraitCor,1)
moduleTraitCor1[abs(moduleTraitCor1)<0.5]=' '

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue1=signif(moduleTraitPvalue,1)
moduleTraitPvalue1[moduleTraitPvalue1<=0.01]='*'
moduleTraitPvalue1[moduleTraitPvalue1>0.01]=' '

sizeGrWindow(20,20)

# genes in different module-------------------
file.remove("All_Gene_KME.txt")
blankdata <- data.frame()
write.table(blankdata,quote = FALSE,row.names = TRUE,file = "All_Gene_KME.txt",col.names = FALSE)
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr0[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote =FALSE,sep='\t',row.names = TRUE,append = TRUE,file = "All_Gene_KME.txt",col.names = FALSE)
}
All_Gene_KME<-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/All_Gene_KME.txt',sep='')
colnames(All_Gene_KME)<-c('id','val','col')
# wgcna.Module-trait.heatmap.pdf-------

MElength=data.frame(length(All_Gene_KME[,1]))
colnames <- colnames(MEs)
colors <- sub("ME", "", colnames)  # 截取"ME"后的部分
#colors <- sub("([A-Za-z]+).*", "\\1", colors)  # 提取颜色部分

# 替换==后的颜色
MElength <- matrix(NA, nrow = length(colors), ncol = 1)

for (i in 1:length(colors)) {
  MElength[i, ] <- length(All_Gene_KME[which(All_Gene_KME$col == colors[i]), 1])
}

  MElength1=as.character(MElength[,1])
  text1 =paste(substring(colnames(MEs[,]),3),"\n",'n=',MElength1)
pdf("wgcna.protein-trait.heatmap.pdf", width =8, height =3)
textMatrix =  paste(moduleTraitCor1[,], "\n",
                   moduleTraitPvalue1[,], "", sep = "");
dim(textMatrix) = dim(moduleTraitCor[,])
par(mar = c(6, 6, 1, 2));
labeledHeatmap(Matrix = moduleTraitCor[,],
               xLabels = names(allTraits),
               yLabels = names(MEs),
               ySymbols = text1,
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
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
#all corlor GO ----

all_plots <- list()
combined_results <- list()
# 循环遍历所有的颜色
for (i in 1:length(colors)) {
  # 选择当前颜色的蛋白数据
  protein_data <- data.frame(All_Gene_KME[which(All_Gene_KME$col == colors[i]), 1])
  proteinlist <- protein_data[, 1]
  
  # 进行GO富集分析
  ego_colors_BP <- enrichGO(gene = proteinlist,
                            keyType = "GID", 
                            OrgDb = Porites.orgdb,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  
  ego_colors_CC <- enrichGO(gene = proteinlist,
                            keyType = "GID", 
                            OrgDb = Porites.orgdb,
                            ont = "CC",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  
  ego_colors_MF <- enrichGO(gene = proteinlist,
                            keyType = "GID", 
                            OrgDb = Porites.orgdb,
                            ont = "MF",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  
  # 提取GO富集结果的数据框
  data_all_colors_BP <- data.frame(ego_colors_BP@result) 
  data_all_colors_CC <- data.frame(ego_colors_CC@result) 
  data_all_colors_MF <- data.frame(ego_colors_MF@result) 
  
  # 添加GO类别列
  data_all_colors_BP$GO <- "BP"
  data_all_colors_CC$GO <- "CC"
  data_all_colors_MF$GO <- "MF"
  
  # 合并结果
  combined_results[[i]] <- rbind(data_all_colors_BP[1:5, ], data_all_colors_CC[1:5, ], data_all_colors_MF[1:5, ])
  
  # 删除重复的Description值
  combined_results[[i]] <- combined_results[[i]][!duplicated(combined_results[[i]]$Description), ]
  
  # 重新创建有序的因子
  combined_results[[i]]$Description <- factor(combined_results[[i]]$Description, levels = unique(combined_results[[i]]$Description), ordered = TRUE)
  
  # 创建并保存图表
  plot <- ggplot(combined_results[[i]])+ ggtitle(colors[i])  +
    aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
    geom_point(shape = 21, color = "black") +
    theme_classic() +
    coord_flip() +
    ylim(0, 30) +
    scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25))) +
    scale_x_discrete(limits = rev(levels(factor(combined_results[[i]]$Description)))) +
    theme(plot.margin = margin(2, 2, 2, 5, "cm"))
  
  all_plots[[i]] <- plot
}

# 保存所有图表到不同的文件
for (i in 1:length(colors)) {
  ggsave(paste("lower7_GO_color_", i, ".pdf", sep = ""), all_plots[[i]], w = 8, h = 10)
}
# different module data-----

# Your original geneid_3 assignment
geneid_3 <- combined_results[[1]]$geneID[11]
# Extract gene names using regular expressions
gene_names <- regmatches(geneid_3, gregexpr("Unigene\\d+", geneid_3))[[1]]

data_imp_norm[gene_names,]
datahgraw_NR[gene_names,]$Nr.annotation


pheatmap(data_imp_norm[gene_names,],scale='row',
         border_color=NA, 
         cluster_cols = FALSE,
         cluster_rows = T,
         col=greenred(75),legend =FALSE)
# 加载 cowplot 包
library(cowplot)

# 创建一个空的图形列表
all_plots_list <- list()

# 循环生成并保存所有图形
for (i in 1:length(colors)) {
  all_plots_list[[i]] <- all_plots[[i]]
}

# 组合所有图形
combined_plot <- plot_grid(plotlist = all_plots_list, ncol = 2, align = 'v')  # 你可以调整 ncol 的值来设置列数

# 保存组合图
ggsave("combined_plots.pdf", combined_plot, w = 20, h = 13)


#all color kegg-------------
# 创建一个空的图形列表
kegg_plots_list <- list()
load("kegg_info.RData")
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/make orgdb/org.P_pukoensis.eg.sqlite")
pathway2gene <- AnnotationDbi::select(Porites.orgdb,
                                      keys = keys(Porites.orgdb),
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

# 循环生成并保存所有颜色的 KEGG 图
for (i in 1:length(colors)){
  # 创建每个颜色的 KEGG 图
  protein_data <- data.frame(All_Gene_KME[which(All_Gene_KME$col == colors[i]), 1])
  proteinlist <- protein_data[, 1]
  
  # 运行 KEGG 富集分析
  ekp_colors_12 <- enricher(proteinlist, TERM2GENE = pathway2gene, TERM2NAME = pathway2name,
                            pvalueCutoff = 0.05, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 5)
  ekp_colors_filtered <- ekp_colors_12@result[ekp_colors_12@result$pvalue < 0.05 & ekp_colors_12@result$Count > 4, ]
  
  # 生成 KEGG 图
  kegg_plot <- ggplot(ekp_colors_filtered[1:10, ]) +
    ggtitle(colors[i]) +
    aes(x = Description, y = -log(pvalue), size = Count) +
    geom_point(shape = 21, color = "black") +
    theme_classic() +
    coord_flip() +
    ylim(0, 20) +
    scale_size_continuous(guide = guide_legend(title = "Count", limits = c(0, 25)))
  
  # 将 KEGG 图添加到列表
  kegg_plots_list[[i]] <- kegg_plot
}

# 组合所有 KEGG 图
combined_kegg_plot <- cowplot::plot_grid(plotlist = kegg_plots_list, ncol = 2, align = 'v')  # 使用 'v' 对齐垂直

# 保存组合的 KEGG 图
ggsave("combined_kegg_plots.pdf", combined_kegg_plot, w = 13, h = 6)
save.image(file = "D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig3WGCNA_GO_KEGG.RData")
