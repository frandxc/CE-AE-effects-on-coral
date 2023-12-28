#BiocManager::install("PerformanceAnalytics")
#Protein and transcript datasets were median centred to the overall median of the respective dataset.

library(ggplot2)
library(ggsignif)
library(vioplot)
library(WGCNA)
library(pheatmap) 
library(RColorBrewer)
library(ggcorrplot)
library(Peptides) 
library("Biostrings")
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')

#all_fasta<- read.fasta('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/Galaxea_Unigene.fa')


load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")

#seqs <- readDNAStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/Galaxea_Unigene.fa', format = "fasta", use.names = TRUE)
fastaFile1<- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/cnidarian/a_digitifera.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
seq_gene=data.frame(pi=pI(fastaFile1))
row.names(seq_gene)=substr(row.names(seq_gene), 1, 14)

aa=as.character(row.names(data_imp_norm) )
p_Nr=datahgraw1[aa,]$Nr.ID
mRNA_Nr=datahgraw1$Nr.ID
seq_p=data.frame(pi=c(seq_gene[p_Nr,]))
seq_mRNA=data.frame(pi=c(seq_gene[mRNA_Nr,]))
#pi density of gene, mRNA, protein------
pdf("Fig4 mRNA_protein.pdf",width=10,height=6)
ggplot() +
  geom_density(data = seq_gene, aes(x = pi, colour = "genome"), alpha = 2, size = 1.5) +
  geom_density(data = seq_mRNA, aes(x = pi, colour = "mRNA"), alpha = 2, size = 1.5)+
  geom_density(data = seq_p, aes(x = pi, colour = "protein"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  scale_x_continuous(breaks=seq(0,15,2))+
  scale_y_continuous(breaks=seq(0,0.23,0.05),limits = c(0,0.23))+
  theme_classic()
dev.off() 


conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
compaired <- list(c("Control", "Ca.OH.2"),
                       c("Control", "CO2"),
                       c("Control", "CO2.Ca.OH.2"),
                       c("Control", "CO2.NaOH"))

row.names(seq_p) <-row.names(data_imp_norm) 

# calculate means
pi_protein<-merge(seq_p,data_imp_norm,by='row.names',all=F)
rownames(pi_protein) <- pi_protein$Row.names
pi_protein_mean <- sapply(seq(1, ncol(pi_protein[,3:22]), by = 4), function(i) {
  rowMeans(pi_protein[,3:22][, i:(i+3)])
})
conditionshg1 <- unique(conditionshg)
colnames(pi_protein_mean)<-conditionshg1 
pi_protein_mean <- data.frame(pi_protein_mean)
pi_protein_mean$Row.names <- rownames(pi_protein_mean) 
pi_protein_mean1 <- merge(pi_protein_mean, pi_protein[, c("Row.names", "pi")], by = "Row.names")
rownames(pi_protein_mean1) <- pi_protein_mean1$Row.names
#pi>7-----
pi_higher_7<-subset(pi_protein_mean1,pi>7)
# 数据矩阵按行进行最小-最大归一化
min_max_normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
normalized_pi_higher_7 <- t(apply(pi_higher_7[,2:6], 1, min_max_normalize))
# 将结果转换为数据框，并设置相同的行名和列名
normalized_pi_higher_7 <- as.data.frame(normalized_pi_higher_7 )
rownames(normalized_pi_higher_7) <- rownames(pi_higher_7[,2:6])
colnames(normalized_pi_higher_7) <- colnames(pi_higher_7[,2:6])

pi_higher_7_to <- tidyr::pivot_longer(normalized_pi_higher_7, cols = colnames(normalized_pi_higher_7), names_to = "Group", values_to = "Y")
group_order <- c("Control", "Ca.OH.2", "CO2", "CO2.Ca.OH.2", "CO2.NaOH")
pi_higher_7_to$Group <- factor(pi_higher_7_to$Group, levels = group_order)

p1= ggplot(pi_higher_7_to, aes(x = Group, y = Y, color = Group)) + 
  labs(y = "relative mean protein abundance") + ggtitle(paste("n =", sum(pi_higher_7_to$Group == "Control"))) +
  theme_classic() +
  coord_cartesian(ylim = c(0.4, 1)) +
  geom_signif(comparisons = compaired,
              map_signif_level = TRUE,
              textsize = 3,
              test = wilcox.test,
              step_increase = 0.1,
              size = 0.5,
              test.args = c("two.sided"),
              y_position = c(0.6, 0.5, 0.4, 0.3)) +
  geom_hline(yintercept = mean(pi_higher_7_to$Y[pi_higher_7_to$Group == "Control"]),
             color = "red", linetype = "dashed", size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black", size = 1) +
  geom_text(aes(x = as.numeric(Group), y = mean(pi_higher_7_to$Y[pi_higher_7_to$Group == "Control"]) + 1,
              label = paste("n =", sum(pi_higher_7_to$Group == "Control"))),
          color = "black")



#pi<7-----

pi_lower_7<-subset(pi_protein_mean1,pi<7)
normalized_pi_lower_7 <- t(apply(pi_lower_7[,2:6], 1, min_max_normalize))
# 将结果转换为数据框，并设置相同的行名和列名
normalized_pi_lower_7 <- as.data.frame(normalized_pi_lower_7 )
rownames(normalized_pi_lower_7) <- rownames(pi_lower_7[,2:6])
colnames(normalized_pi_lower_7) <- colnames(pi_lower_7[,2:6])
pi_lower_7_to <- tidyr::pivot_longer(normalized_pi_lower_7, cols = colnames(normalized_pi_lower_7), names_to = "Group", values_to = "Y")
group_order <- c("Control", "Ca.OH.2", "CO2", "CO2.Ca.OH.2", "CO2.NaOH")
pi_lower_7_to$Group <- factor(pi_lower_7_to$Group, levels = group_order)

# 绘制带点的箱线图
p2= ggplot(pi_lower_7_to, aes(x = Group, y = Y, color = Group)) +
  labs(y = "relative mean protein abundance") + ggtitle(paste("n =", sum(pi_lower_7_to$Group == "Control"))) +
  theme_classic() +
  coord_cartesian(ylim = c(0.4, 1)) +
  geom_signif(comparisons = compaired,
              map_signif_level = TRUE,
              textsize = 3,
              test = wilcox.test,
              step_increase = 0.1,
              size = 0.5,
              test.args = c("two.sided"),
              y_position = c(0.6, 0.5, 0.4, 0.3)) +
  geom_hline(yintercept = mean(pi_lower_7_to$Y[pi_lower_7_to$Group == "Control"]),
             color = "red", linetype = "dashed", size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black", size = 1) +
  geom_text(aes(x = as.numeric(Group), y = mean(pi_lower_7_to$Y[pi_lower_7_to$Group == "Control"]) + 1,
                label = paste("n =", sum(pi_lower_7_to$Group == "Control"))),
            color = "black")
library(gridExtra)

# 创建示例的两个ggplot对象 p1 和 p2

# 将p1和p2按照两列合并


# 保存合并后的图形到PDF文件
pdf("combined_plots.pdf", width = 10, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

datahgraw_NR1 <-datahgraw1
# carbonic -----
carbonic <- datahgraw_NR[grep(pattern = "Integrin", datahgraw_NR$Nr.annotation), ]
p_carbonic_heatmap<-normalized_hg[rownames(normalized_hg) %in% rownames(carbonic), ]
p_carbonic_heatmap1<-p_carbonic_heatmap[rowSums(p_carbonic_heatmap[,1:20])>1,]
p_carbonic_heatmap_g<-data.frame(p_carbonic_heatmap[rownames(p_carbonic_heatmap) %in% rownames(p_carbonic_heatmap1), ])
p_carbonic_heatmap_g[,21:(20+ncol(datahgraw_NR))]<-datahgraw_NR[row.names(p_carbonic_heatmap_g),]

p_carbonic_heatmap_g_labels <- p_carbonic_heatmap_g$Nr.annotation
max_label_length <-30
# 将标签分割成固定长度的片段，并创建包含换行符的字符串
p_carbonic_heatmap_g_labels <- sapply(p_carbonic_heatmap_g_labels, function(label) {
  paste(strwrap(label, width = max_label_length), collapse = "\n")
})
pheatmap(
  p_carbonic_heatmap_g[,1:20],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
  labels_row = p_carbonic_heatmap_g_labels,
  labels_col = conditionshg,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)

#pi ---
carbonic1 <- datahgraw_NR[rownames(carbonic),]$Nr.ID
seq_gene$ID  <- rownames(seq_gene)
carbonic2 <-seq_gene[rownames(seq_gene) %in% carbonic1,]




seq_gene=data.frame(pi=pI(fastaFile1))


aa=as.character(row.names(data_imp_norm) )
p_Nr=datahgraw1[aa,]$Nr.ID
mRNA_Nr=datahgraw1$Nr.ID
seq_p=data.frame(pi=c(seq_gene[p_Nr,]))
seq_mRNA=data.frame(pi=c(seq_gene[mRNA_Nr,]))