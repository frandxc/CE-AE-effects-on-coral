#FIG  protein enrichment of 1) pi >7 <7 ; 2) correlated with CO2;2) correlated with alkalinity
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
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/make orgdb/org.P_pukoensis.eg.sqlite")
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig.4_PI7.RData")
#GO BP >7 -----
protein_data_higher7=subset(pi_protein,pi>7)
ego_pi_higher7_BP <- enrichGO(gene = protein_data_higher7[,1],
                              keyType = "GID", 
                              OrgDb =Porites.orgdb,
                              ont = "BP", #ALL BP或MF或CC
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05) 
ego_pi_higher7_BP_f <-ego_pi_higher7_BP@result[ego_pi_higher7_BP@result$pvalue <0.05&ego_pi_higher7_BP@result$Count>3,]

#GO CC higher7-----
ego_pi_higher7_CC <- enrichGO(gene = protein_data_higher7[,1],
                              keyType = "GID", 
                              OrgDb =Porites.orgdb,
                              ont = "CC", #ALL BP或MF或CC
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05) 
ego_pi_higher7_CC_f <-ego_pi_higher7_CC@result[ego_pi_higher7_CC@result$pvalue <0.05&ego_pi_higher7_CC@result$Count>3,]
#GO CC higher7-----
ego_pi_higher7_CC <- enrichGO(gene = protein_data_higher7[,1],
                              keyType = "GID", 
                              OrgDb =Porites.orgdb,
                              ont = "CC", #ALL BP或MF或CC
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05) 
ego_pi_higher7_CC_f <-ego_pi_higher7_CC@result[ego_pi_higher7_CC@result$pvalue <0.05&ego_pi_higher7_CC@result$Count>3,]
#GO MF higher7-----
ego_pi_higher7_MF <- enrichGO(protein_data_higher7[,1], 
                              keyType = "GID", 
                              OrgDb =Porites.orgdb,
                              ont = "MF", #ALL BP或MF或CC
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05) 
ego_pi_higher7_MF_f <- ego_pi_higher7_MF@result[ego_pi_higher7_MF@result$pvalue <0.05&ego_pi_higher7_MF@result$Count>3,]
ego_pi_higher7_BP_f$GO <- "BP"
ego_pi_higher7_CC_f$GO <- "CC"
ego_pi_higher7_MF_f$GO <- "MF"
pdata_higher7 = rbind(ego_pi_higher7_BP_f[1:5,],ego_pi_higher7_CC_f[1:5,],ego_pi_higher7_MF_f[1:5,])
pdata_higher7$Description<- factor(pdata_higher7$Description, level=pdata_higher7$Description, ordered=TRUE)

higher7_GO <- ggplot(pdata_higher7) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape = 21, color = "black") +
  theme_classic() +
  coord_flip() +
  ylim(0, 20) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25))) +
  scale_x_discrete(limits = rev(levels(factor(pdata_higher7$Description)))) +
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))

ggsave("higher7_GO.pdf", higher7_GO,  w=8, h=10)


# 
# dotplot(ego_pi_7.5_ALL ,x = "Count", color = "pvalue", size = "GeneRatio", #默认参数
#         split="ONTOLOGY") + #以ONTOLOGY类型分开
#   facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分屏绘图


#GO BP <7-------
protein_data_7=subset(pi_protein,pi<7)

ego_pi_7_BP <- enrichGO(gene = protein_data_7[,1],
                        keyType = "GID", 
                        OrgDb =Porites.orgdb,
                        ont = "BP", #ALL BP或MF或CC
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05) 
ego_pi_7_BP_f <-ego_pi_7_BP@result[ego_pi_7_BP@result$pvalue <0.05&ego_pi_7_BP@result$Count>3,]
#GO CC <7-------
ego_pi_7_CC <- enrichGO(gene =protein_data_7[,1],
                        keyType = "GID", 
                        OrgDb =Porites.orgdb,
                        ont = "CC", #ALL BP或MF或CC
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05) 
ego_pi_7_CC_f <-ego_pi_7_CC@result[ego_pi_7_CC@result$pvalue <0.05&ego_pi_7_CC@result$Count>3,]
#GO MF <7-------
ego_pi_7_MF <- enrichGO(gene =protein_data_7[,1],
                        keyType = "GID", 
                        OrgDb =Porites.orgdb,
                        ont = "MF", #ALL BP或MF或CC
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05) 
ego_pi_7_MF_f <-ego_pi_7_MF@result[ego_pi_7_MF@result$pvalue <0.05&ego_pi_7_MF@result$Count>3,]

ego_pi_7_BP_f$GO <- "BP"
ego_pi_7_CC_f$GO <- "CC"
ego_pi_7_MF_f$GO <- "MF"
pdata_lower7 = rbind(ego_pi_7_BP_f[1:5,],ego_pi_7_CC_f[1:5,],ego_pi_7_MF_f[1:5,])
pdata_lower7$Description<- factor(pdata_lower7$Description, level=pdata_lower7$Description, ordered=TRUE)

lower7_GO <-ggplot(pdata_lower7) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape = 21, color = "black") +
  theme_classic() +
  coord_flip() +
  ylim(0, 20) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25))) +
  scale_x_discrete(limits = rev(levels(factor(pdata_lower7$Description)))) +
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))
ggsave("lower7_GO.pdf", lower7_GO,  w=8, h=10)
#all plot -----
library(cowplot)
# 组合所有图形并设置主图的尺寸
pdf("all_go_highandlow.pdf", width =20, height =8)
plot_grid(lower7_GO, higher7_GO,
          ncol = 2, align = "v",
          labels = c("a", "b",
                     label_size = 14,
                     rel_heights = c(1, 1),
                     rel_widths = c(1, 1)))
dev.off()

# positive AND NEG correlate with CO2------
traits <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/traits.csv', header = T, row.names=1) 
ALL_cor<- data.frame(condition=conditionshg,t(data_imp_norm[,1:20]),CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)
cor_all<- data.frame(cor(ALL_cor[,-1]))

cor_all1<-cor_all[, (length(cor_all[,1]) - 4):length(cor_all[,1])]
cor_all_p <- cor_all1[cor_all1$CO2>0.5,]
cor_all_n <- cor_all1[cor_all1$CO2 < -0.5,]

Porites.orgdb <- loadDb("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/make orgdb/org.P_pukoensis.eg.sqlite")
#positive WITH CO2----
proteinlist2 = rownames(cor_all_p)
ego_positive_co2_BP <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "BP", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_co2_GO <-data.frame(ego_positive_co2_BP@result) 

ego_positive_co2_CC <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "CC", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_co2_CC <-data.frame(ego_positive_co2_CC@result) 
ego_positive_co2_MF <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "MF", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_co2_MF <-data.frame(ego_positive_co2_MF@result) 
data_all_positive_co2_MF$GO <- "MF"
data_all_positive_co2_CC$GO <- "CC"
data_all_positive_co2_GO$GO <- "BP"
pdata = rbind(data_all_positive_co2_GO[1:5,],data_all_positive_co2_CC[1:5,],data_all_positive_co2_MF[1:5,])
pdata$Description<- factor(pdata$Description, level=pdata$Description, ordered=TRUE)

head(pdata)
positve_GO <- ggplot(pdata) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape = 21, color = "black")+ theme_bw() +coord_flip() +ylim(0,20) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25)))+
  scale_x_discrete(limits = rev(levels(factor(pdata$Description))))+
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))



# negative WITH CO2----
proteinlist_neg = rownames(cor_all_n)
ego_negative_BP <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "BP", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_GO <-data.frame(ego_negative_BP@result) 

ego_negative_CC <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "CC", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_CC <-data.frame(ego_negative_CC@result) 
ego_negative_MF <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "MF", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_MF <-data.frame(ego_negative_MF@result) 
data_all_negative_MF$GO <- "MF"
data_all_negative_CC$GO <- "CC"
data_all_negative_GO$GO <- "BP"
pdata_negative = rbind(data_all_negative_GO[1:5,],data_all_negative_CC[1:5,],data_all_negative_MF[1:5,])
pdata_negative$Description<- factor(pdata_negative$Description, level=pdata_negative$Description, ordered=TRUE)

negative_GO <- ggplot(pdata_negative) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape =21, color = "black")+ theme_bw() +coord_flip()+ylim(0,20)+ 
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25)))+
  scale_x_discrete(limits = rev(levels(factor(pdata_negative$Description))))+
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))

#all plot -----
library(cowplot)
# 组合所有图形并设置主图的尺寸
pdf("all_go.pdf", width =20, height =16)
plot_grid(positve_GO,negative_GO, phigher7,plower7,
          ncol = 2, align = "v",
          labels = c("A)", "B)", "C)", "D)"),
          label_size = 14,
          rel_heights = c(1, 1),
          rel_widths = c(1, 1))
dev.off()
# positive AND NEG correlate with alk------
cor_alk_p <- cor_all1[cor_all1$alk>0.5,]
cor_alk_n <- cor_all1[cor_all1$alk < -0.5,]

#positive WITH alk----
proteinlist2 = rownames(cor_alk_p)
ego_positive_alk_BP <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "BP", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_alk_GO <-data.frame(ego_positive_alk_BP@result) 

ego_positive_alk_CC <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "CC", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_alk_CC <-data.frame(ego_positive_alk_CC@result) 
ego_positive_alk_MF <- enrichGO(gene = proteinlist2,
                                keyType = "GID", 
                                OrgDb =Porites.orgdb,
                                ont = "MF", #ALL BP或MF或CC
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05) 
data_all_positive_alk_MF <-data.frame(ego_positive_alk_MF@result) 
data_all_positive_alk_MF$GO <- "MF"
data_all_positive_alk_CC$GO <- "CC"
data_all_positive_alk_GO$GO <- "BP"
pdata = rbind(data_all_positive_alk_GO[1:5,],data_all_positive_alk_CC[1:5,],data_all_positive_alk_MF[1:5,])
pdata$Description<- factor(pdata$Description, level=pdata$Description, ordered=TRUE)

head(pdata)
positve_GO_alk <- ggplot(pdata) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape = 21, color = "black")+ theme_bw() +coord_flip() +ylim(0,20) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25)))+
  scale_x_discrete(limits = rev(levels(factor(pdata$Description))))+
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))



# negative WITH alk----
proteinlist_neg = rownames(cor_alk_n)
ego_negative_BP <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "BP", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_GO <-data.frame(ego_negative_BP@result) 

ego_negative_CC <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "CC", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_CC <-data.frame(ego_negative_CC@result) 
ego_negative_MF <- enrichGO(gene = proteinlist_neg,
                            keyType = "GID", 
                            OrgDb =Porites.orgdb,
                            ont = "MF", #ALL BP或MF或CC
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05) 
data_all_negative_MF <-data.frame(ego_negative_MF@result) 
data_all_negative_MF$GO <- "MF"
data_all_negative_CC$GO <- "CC"
data_all_negative_GO$GO <- "BP"
pdata_negative = rbind(data_all_negative_GO[1:5,],data_all_negative_CC[1:5,],data_all_negative_MF[1:5,])
pdata_negative$Description<- factor(pdata_negative$Description, level=pdata_negative$Description, ordered=TRUE)

negative_GO_alk <- ggplot(pdata_negative) +
  aes(x = Description, y = -log(pvalue), fill = GO, size = Count) +
  geom_point(shape =21, color = "black")+ theme_bw() +coord_flip()+ylim(0,20)+ 
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 25)))+
  scale_x_discrete(limits = rev(levels(factor(pdata_negative$Description))))+
  theme(plot.margin = margin(2, 2, 2, 5, "cm"))

#all plot -----
library(cowplot)
# 组合所有图形并设置主图的尺寸
pdf("all_go_alk.pdf", width =20, height =8)
plot_grid(positve_GO_alk,negative_GO_alk,
          ncol = 2, align = "v",
          labels = c("A)", "B)",
                     label_size = 14,
                     rel_heights = c(1, 1),
                     rel_widths = c(1, 1)))
dev.off()

# GO subcellular--------
# mitocho-----
seq_gene$Nr.ID =rownames(seq_gene)
mitocho <- datahgraw1[grep(pattern = "mitochondrial", datahgraw1$GO.CellularComponent), ]
mitocho_m <- merge(seq_gene,mitocho,by="Nr.ID")

cytoplasmic <- datahgraw1[grep(pattern = "cytoplasm", datahgraw1$GO.CellularComponent), ]
cytoplasmic <- merge(seq_gene,cytoplasmic,by="Nr.ID")

ribosome <- datahgraw1[grep(pattern = "ribosome", datahgraw1$GO.CellularComponent), ]
ribosome_m <- merge(seq_gene,ribosome,by="Nr.ID")

nucl <- datahgraw1[grep(pattern = "nucleus", datahgraw1$GO.CellularComponent), ]
nucl_m <- merge(seq_gene,nucl,by="Nr.ID")

reticulum <- datahgraw1[grep(pattern = "reticulum", datahgraw1$GO.CellularComponent), ]
reticulum_m <- merge(seq_gene,reticulum,by="Nr.ID")

Golgi <- datahgraw1[grep(pattern = "Golgi apparatus", datahgraw1$GO.CellularComponent), ]
Golgi_m <- merge(seq_gene,Golgi,by="Nr.ID")


membrane <- datahgraw1[grep(pattern = "membrane", datahgraw1$GO.CellularComponent), ]
membrane_m <- merge(seq_gene,membrane,by="Nr.ID")


organelle_inner  <- datahgraw1[grep(pattern = "organelle inner membrane", datahgraw1$GO.CellularComponent), ]
organelle_inner_m <- merge(seq_gene,organelle_inner,by="Nr.ID")

endosome  <- datahgraw1[grep(pattern = "endosome", datahgraw1$GO.CellularComponent), ]
endosome_m <- merge(seq_gene,endosome,by="Nr.ID")

peroxisome   <- datahgraw1[grep(pattern = "peroxisome", datahgraw1$GO.CellularComponent), ]
peroxisome_m <- merge(seq_gene,peroxisome ,by="Nr.ID")

# 创建一个空的数据框用于存储结果
result_df <- data.frame(Group = character(),
                        Mean = numeric(),
                        SD = numeric(),
                        stringsAsFactors = FALSE)

# 定义一个函数用于计算平均值和方差，并将结果添加到结果数据框中
calculate_mean_sd <- function(group_name, group_data) {
  mean_value <- mean(group_data$pi)
  sd_value <- sd(group_data$pi)
  
  result_df <<- rbind(result_df, data.frame(Group = group_name,
                                            Mean = mean_value,
                                            SD = sd_value,
                                            stringsAsFactors = FALSE))
}

# 分别计算每个细胞组分的平均值和方差，并将结果添加到结果数据框中
calculate_mean_sd("mitochondrial", mitocho_m)
calculate_mean_sd("cytoplasmic", cytoplasmic)
calculate_mean_sd("ribosome", ribosome_m)
calculate_mean_sd("nucleus", nucl_m)
calculate_mean_sd("endosome", endosome_m)
calculate_mean_sd("Golgi apparatus", Golgi_m)
calculate_mean_sd("organelle inner membrane", organelle_inner_m)
calculate_mean_sd("reticulum", reticulum_m)
calculate_mean_sd("membrane", membrane_m)
calculate_mean_sd("peroxisome", peroxisome_m)

# 查看结果数据框
result_df

library(ggplot2)

# 绘制平均值的柱状图
mean_plot <- ggplot(result_df, aes(x = Group, y = Mean)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.5) +
  labs(x = "Cellular Component", y = "Mean pI") +
  theme_minimal()

# 绘制方差的误差棒图
sd_plot <- ggplot(result_df, aes(x = Group, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_errorbar(width = 0.2, color = "black") +
  labs(x = "Cellular Component", y = "Mean pI") +
  theme_minimal()

# 合并两个图层
combined_plot <- mean_plot + sd_plot

# 显示图形
combined_plot

# 选择与cytoplasmic组对比的数据
comparison_df <- result_df[result_df$Group != "cytoplasmic", ]

# 绘制误差棒图和散点图
subcellular <- ggplot(comparison_df, aes(x = Group, y = Mean, color = Group)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = '', y = 'pI')+  geom_hline(yintercept = 7,
                                      color = "red", linetype = "dashed", size = 1)

ggsave("subcellular.pdf",subcellular , w=6, h=6)
