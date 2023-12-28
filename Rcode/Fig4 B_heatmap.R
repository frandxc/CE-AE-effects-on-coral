library(ggplot2)
library("Biostrings")
library(Peptides) 
library(rentrez)
library(pheatmap)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")

#heat_CaOH2_up----
heat_CaOH2_up<-data_imp_coral[data_imp_coral$Protein %in% CaOH2_up$name, ][,c(1,4:23)]

P11<-pheatmap(
  heat_CaOH2_up[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CaOH2_down----
heat_CaOH2_down<-data_imp_coral[data_imp_coral$Protein %in% CaOH2_down$name, ][,c(1,4:23)]
P12<-pheatmap(
  heat_CaOH2_down[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2_up----
heat_CO2_up<-data_imp_coral[data_imp_coral$Protein %in% CO2_up$name, ][,c(1,4:23)]
P13<-pheatmap(
  heat_CO2_up[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2_down----
heat_CO2_down<-data_imp_coral[data_imp_coral$Protein %in% CO2_down$name, ][,c(1,4:23)]
P14<-pheatmap(
  heat_CO2_down[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2CaOH2_up----
heat_CO2CaOH2_up<-data_imp_coral[data_imp_coral$Protein %in% CO2CaOH2_up$name, ][,c(1,4:23)]
P15<-pheatmap(
  heat_CO2CaOH2_up[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2CaOH2_down----
heat_CO2CaOH2_down<-data_imp_coral[data_imp_coral$Protein %in% CO2CaOH2_down$name, ][,c(1,4:23)]
P16<-pheatmap(
  heat_CO2CaOH2_down[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2NaOH_up----
heat_CO2NaOH_up<-data_imp_coral[data_imp_coral$Protein %in% CO2NaOH_up$name, ][,c(1,4:23)]
P17<-pheatmap(
  heat_CO2NaOH_up[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)
#heat_CO2NaOH_down----
heat_CO2NaOH_down<-data_imp_coral[data_imp_coral$Protein %in% CO2NaOH_down$name, ][,c(1,4:23)]
P18<-pheatmap(
  heat_CO2NaOH_down[,2:21],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = F,
  legend = FALSE,
    labels_row = "",
  labels_col = "",
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)

#  fasta------
library(Peptides) #calculate pI
library(rentrez) #extract sequence
# # 定义处理的条件
# conditions <- c("CaOH2_up", "CaOH2_down", "CO2_up", "CO2_down", "CO2CaOH2_up", "CO2CaOH2_down","CO2NaOH_up", "CO2NaOH_down")
# 
# # 循环遍历每个条件
# for (condition in conditions) {
#   # 获取相应条件的Accession号
#   condition_name <- get(condition)
#   accession_numbers <- datahgraw[datahgraw$geneID %in% condition_name$name, ]$Nr.ID[datahgraw[datahgraw$geneID %in% condition_name$name, ]$Nr.ID != "--"]
#   
#   # 初始化字符向量
#   protein_sequences <- character(0)
#   
#   # 循环遪Accession号并下载蛋白序列
#   for (accession in accession_numbers) {
#     # 使用efetch函数下载蛋白序列
#     protein_seq <- entrez_fetch(db = "protein", id = accession, rettype = "fasta", retmode = "text")
#     # 打印蛋白序列
#     cat(protein_seq, "\n")
#     # 将蛋白序列附加到字符向量中
#     protein_sequences <- c(protein_sequences, protein_seq)
#   }
#   
#   # 保存蛋白序列到相应的文件
#   write(paste(protein_sequences, collapse = "\n"), file = paste0(condition, "_protein_sequences.fasta"))
# }

fastaFile_CaOH2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CaOH2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CaOH2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CaOH2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2CaOH2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2CaOH2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2CaOH2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2CaOH2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2NaOH_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2NaOH_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2NaOH_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2NaOH_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)

#CaOH2 UP------
fastaFile_CaOH2_up_name1 <- names(fastaFile_CaOH2_up) %>% as.data.frame()
fastaFile_CaOH2_up_data1 <- as.data.frame(fastaFile_CaOH2_up)
fastaFile_CaOH2_up_data21 <- data.frame(ID= fastaFile_CaOH2_up_name1$.,seq=fastaFile_CaOH2_up_data1$x, ID2= substr( fastaFile_CaOH2_up_name1$., 1, 14))
fastaFile_CaOH2_up_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CaOH2_up_data21$ID2)
fastaFile_CaOH2_up_seq1=data.frame(fastaFile_CaOH2_up_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CaOH2_up_data21$seq),1)) 
fastaFile_CaOH2_up_seq2=unique(fastaFile_CaOH2_up_seq1, by = "ID3")
P1<-ggplot(fastaFile_CaOH2_up_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CaOH2_up_seq2, PI1 < 7))
nrow(subset(fastaFile_CaOH2_up_seq2, PI1 > 7))

#CaOH2 down------
fastaFile_CaOH2_down_name1 <- names(fastaFile_CaOH2_down) %>% as.data.frame()
fastaFile_CaOH2_down_data1 <- as.data.frame(fastaFile_CaOH2_down)
fastaFile_CaOH2_down_data21 <- data.frame(ID= fastaFile_CaOH2_down_name1$.,seq=fastaFile_CaOH2_down_data1$x, ID2= substr( fastaFile_CaOH2_down_name1$., 1, 14))
fastaFile_CaOH2_down_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CaOH2_down_data21$ID2)
fastaFile_CaOH2_down_seq1=data.frame(fastaFile_CaOH2_down_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CaOH2_down_data21$seq),1)) 
fastaFile_CaOH2_down_seq2=unique(fastaFile_CaOH2_down_seq1, by = "ID3")
P2<-ggplot(fastaFile_CaOH2_down_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CaOH2_down_seq2, PI1 < 7))
nrow(subset(fastaFile_CaOH2_down_seq2, PI1 > 7))
#CO2 UP------
fastaFile_CO2_up_name1 <- names(fastaFile_CO2_up) %>% as.data.frame()
fastaFile_CO2_up_data1 <- as.data.frame(fastaFile_CO2_up)
fastaFile_CO2_up_data21 <- data.frame(ID= fastaFile_CO2_up_name1$.,seq=fastaFile_CO2_up_data1$x, ID2= substr( fastaFile_CO2_up_name1$., 1, 14))
fastaFile_CO2_up_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2_up_data21$ID2)
fastaFile_CO2_up_seq1=data.frame(fastaFile_CO2_up_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2_up_data21$seq),1)) 
fastaFile_CO2_up_seq2=unique(fastaFile_CO2_up_seq1, by = "ID3")
P3<-ggplot(fastaFile_CO2_up_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2_up_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2_up_seq2, PI1 > 7))
#CO2 down------
fastaFile_CO2_down_name1 <- names(fastaFile_CO2_down) %>% as.data.frame()
fastaFile_CO2_down_data1 <- as.data.frame(fastaFile_CO2_down)
fastaFile_CO2_down_data21 <- data.frame(ID= fastaFile_CO2_down_name1$.,seq=fastaFile_CO2_down_data1$x, ID2= substr( fastaFile_CO2_down_name1$., 1, 14))
fastaFile_CO2_down_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2_down_data21$ID2)
fastaFile_CO2_down_seq1=data.frame(fastaFile_CO2_down_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2_down_data21$seq),1)) 
fastaFile_CO2_down_seq2=unique(fastaFile_CO2_down_seq1, by = "ID3")
P4<-ggplot(fastaFile_CO2_down_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2_down_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2_down_seq2, PI1 > 7))
#CO2CaOH2 UP------
fastaFile_CO2CaOH2_up_name1 <- names(fastaFile_CO2CaOH2_up) %>% as.data.frame()
fastaFile_CO2CaOH2_up_data1 <- as.data.frame(fastaFile_CO2CaOH2_up)
fastaFile_CO2CaOH2_up_data21 <- data.frame(ID= fastaFile_CO2CaOH2_up_name1$.,seq=fastaFile_CO2CaOH2_up_data1$x, ID2= substr( fastaFile_CO2CaOH2_up_name1$., 1, 14))
fastaFile_CO2CaOH2_up_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2CaOH2_up_data21$ID2)
fastaFile_CO2CaOH2_up_seq1=data.frame(fastaFile_CO2CaOH2_up_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2CaOH2_up_data21$seq),1)) 
fastaFile_CO2CaOH2_up_seq2=unique(fastaFile_CO2CaOH2_up_seq1, by = "ID3")
P5<-ggplot(fastaFile_CO2CaOH2_up_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2CaOH2_up_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2CaOH2_up_seq2, PI1 > 7))
#CO2CaOH2 down------
fastaFile_CO2CaOH2_down_name1 <- names(fastaFile_CO2CaOH2_down) %>% as.data.frame()
fastaFile_CO2CaOH2_down_data1 <- as.data.frame(fastaFile_CO2CaOH2_down)
fastaFile_CO2CaOH2_down_data21 <- data.frame(ID= fastaFile_CO2CaOH2_down_name1$.,seq=fastaFile_CO2CaOH2_down_data1$x, ID2= substr( fastaFile_CO2CaOH2_down_name1$., 1, 14))
fastaFile_CO2CaOH2_down_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2CaOH2_down_data21$ID2)
fastaFile_CO2CaOH2_down_seq1=data.frame(fastaFile_CO2CaOH2_down_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2CaOH2_down_data21$seq),1)) 
fastaFile_CO2CaOH2_down_seq2=unique(fastaFile_CO2CaOH2_down_seq1, by = "ID3")
P6<-ggplot(fastaFile_CO2CaOH2_down_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2CaOH2_down_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2CaOH2_down_seq2, PI1 > 7))
#CO2NaOH UP------
fastaFile_CO2NaOH_up_name1 <- names(fastaFile_CO2NaOH_up) %>% as.data.frame()
fastaFile_CO2NaOH_up_data1 <- as.data.frame(fastaFile_CO2NaOH_up)
fastaFile_CO2NaOH_up_data21 <- data.frame(ID= fastaFile_CO2NaOH_up_name1$.,seq=fastaFile_CO2NaOH_up_data1$x, ID2= substr( fastaFile_CO2NaOH_up_name1$., 1, 14))
fastaFile_CO2NaOH_up_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2NaOH_up_data21$ID2)
fastaFile_CO2NaOH_up_seq1=data.frame(fastaFile_CO2NaOH_up_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2NaOH_up_data21$seq),1)) 
fastaFile_CO2NaOH_up_seq2=unique(fastaFile_CO2NaOH_up_seq1, by = "ID3")
P7<-ggplot(fastaFile_CO2NaOH_up_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic()  +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2NaOH_up_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2NaOH_up_seq2, PI1 > 7))
#CO2NaOH down------
fastaFile_CO2NaOH_down_name1 <- names(fastaFile_CO2NaOH_down) %>% as.data.frame()
fastaFile_CO2NaOH_down_data1 <- as.data.frame(fastaFile_CO2NaOH_down)
fastaFile_CO2NaOH_down_data21 <- data.frame(ID= fastaFile_CO2NaOH_down_name1$.,seq=fastaFile_CO2NaOH_down_data1$x, ID2= substr( fastaFile_CO2NaOH_down_name1$., 1, 14))
fastaFile_CO2NaOH_down_data21$ID3<-  sub("(\\.\\d).*", "\\1", fastaFile_CO2NaOH_down_data21$ID2)
fastaFile_CO2NaOH_down_seq1=data.frame(fastaFile_CO2NaOH_down_data21) %>%
  mutate("PI1" = round(pI(seq = fastaFile_CO2NaOH_down_data21$seq),1)) 
fastaFile_CO2NaOH_down_seq2=unique(fastaFile_CO2NaOH_down_seq1, by = "ID3")
P8<-ggplot(fastaFile_CO2NaOH_down_seq2, aes(x = PI1, y =  reorder(ID3, -PI1) )) +
  geom_bar(stat = "identity", width = 0.5, fill = "blue") +
  theme_classic() +  xlim(0,15)+
  geom_vline(xintercept = 7, linetype = "dashed", color = "red")+ theme(axis.text.y = element_blank())    # 在x=7处画一条虚线
nrow(subset(fastaFile_CO2NaOH_down_seq2, PI1 < 7))
nrow(subset(fastaFile_CO2NaOH_down_seq2, PI1 > 7))
heat_plots <- plot_grid(P11$gtable,P12$gtable,P13$gtable,P14$gtable,P15$gtable,P16$gtable,P17$gtable,P18$gtable, ncol = 2, align = "hv")  # "hv" stands for horizontal  
pi_bar_plots <- plot_grid(P1,P2,P3,P4,P5,P6,P7,P8, ncol = 2, align = "hv")  # "hv" stands for horizontal  
ggsave("pI_combined_heatmaps.pdf", heat_plots, width =8, height = 12)
ggsave("pI_combined_bar.pdf", pi_bar_plots, width = 8, height = 12)
