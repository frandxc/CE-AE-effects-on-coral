setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data')
library(VennDiagram)
library(UpSetR)
library(ggvenn)
library(forcats)
require(cowplot)
library(gridExtra)
library(ggplot2)
library(stringr)
library(Peptides)  
library(tidyverse)
library(pheatmap) 
library(rentrez)
library(stats)
library(ggsignif)

genome_ex <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/adi_hydra_desmo/expanded_GO.csv', header = F) 
genome_con <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/adi_hydra_desmo/contracted_GO.csv', header = F)
# split columns of genome_ex
split_columns <- strsplit(as.character(genome_ex$V1), "\\s+")
genome_ex$GOID <- sapply(split_columns, function(x) x[1])
genome_ex$count <- sapply(split_columns, function(x) as.numeric(x[2]))
genome_ex$description <- sapply(split_columns, function(x) paste(x[3:(length(x)-2)], collapse = " "))
genome_ex$GO_category <- sapply(split_columns, function(x) x[length(x)-1])
genome_ex$evalue <- sapply(split_columns, function(x) as.numeric(x[length(x)]))
genome_ex$V1 <- NULL
# split columns of genome_con
split_columns1 <- strsplit(as.character(genome_con$V1), "\\s+")
genome_con$GOID <- sapply(split_columns1, function(x) x[1])
genome_con$count <- sapply(split_columns1, function(x) as.numeric(x[2]))
genome_con$description <- sapply(split_columns1, function(x) paste(x[3:(length(x)-2)], collapse = " "))
genome_con$GO_category <- sapply(split_columns1, function(x) x[length(x)-1])
genome_con$evalue <- sapply(split_columns1, function(x) as.numeric(x[length(x)]))
genome_con$V1 <- NULL


# fig expanded------
p1= genome_ex[1:10,] %>%
  mutate(name = fct_reorder(description, -log(evalue))) %>%
  ggplot(aes(x = name, y = -log(evalue), fill = GO_category, size =count) )+
  geom_point(shape = 21, color = "black") +
  coord_flip() +
  ylim(0, 100) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 30))) +
  theme_classic()
# fig contracted------
genome_con1 <- genome_con %>%
  mutate(name = fct_reorder(description, count))
p2= genome_con1[1:10,] %>%
  ggplot(aes(x = name, y = -log(evalue), fill = GO_category, size =count) )+
  geom_point(shape = 21, color = "black") +
  coord_flip() +
  ylim(0, 100) +
  scale_size_continuous(range = c(1, 15), guide = guide_legend(title = "Count", limits = c(0, 30))) +
  theme_classic()

library(patchwork)
ex_con <- p1/p2+plot_layout(widths = c(1, 1))
ggsave("ex_con.pdf", ex_con,  w=12, h=10)


# Cnidarian PI and sequence------
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/cnidarian')
datahgraw1 <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/All_Unigene.advance.annotation.csv', header = T, row.names=1) 
ex_raw1 <-read.csv("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/adi_hydra_desmo/expanded_cluster.csv",  header = F) #expanded gene family
head(ex_raw1)
con_raw1 <- read.csv("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/adi_hydra_desmo/contracted_cluster.csv",  header = F)#contracted gene family

a <- list.files("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/cnidarian")  
a
fastaFile1<- Biostrings::readBStringSet('a_digitifera.fasta', format = "fasta", use.names = TRUE)
fastaFile2<- Biostrings::readBStringSet('Hydra vulgaris fresh water hydra.fasta', format = "fasta", use.names = TRUE)
fastaFile3<- Biostrings::readBStringSet('Desmophyllum_pertusum_deepsea.faa', format = "fasta", use.names = TRUE)
# digitifera_data <- read.table("a_digitifera.csv", sep = "\t")
# Hydra_data <- read.table("Hydra_vulgaris.csv", sep = "\t")
# 使用substr()函数拆分列并创建两个新列
# digitifera_data$ID2 <- substr(digitifera_data$V1, 1, 14)
# digitifera_data$ANN <- substr(digitifera_data$V1, 16, nchar(digitifera_data$V1,))
# digitifera_data$V1 <- NULL
# Hydra_data$ID2 <-substr(Hydra_data$V1, 1, 14)
# Hydra_data$ANN <-substr(Hydra_data$V1, 16, nchar(Hydra_data$V1,))
# Hydra_data$V1 <-NULL
# 打印数据框内容
name1 <- names(fastaFile1) %>% as.data.frame()
name2 <- names(fastaFile2) %>% as.data.frame()
name3 <- names(fastaFile3) %>% as.data.frame()

data1 <- as.data.frame(fastaFile1)
data2 <- as.data.frame(fastaFile2)
data3 <- as.data.frame(fastaFile3)

data21 <- data.frame(ID= name1$.,seq=data1$x, ID2= substr( name1$., 1, 14))
data22 <- data.frame(ID= name2$.,seq=data2$x, ID2= substr( name2$., 1, 14))
data23 <- data.frame(ID= name3$.,seq=data3$x, ID2= substr( name3$., 1, 12))

seq1=data.frame(data21) %>%
  mutate("PI1" = pI(seq = data21$seq)) 
seq2=data.frame(data22) %>%
  mutate("PI1" = pI(seq = data22$seq)) 
seq3=data.frame(data23) %>%
  mutate("PI1" = pI(seq = data23$seq)) 
all_seq <- rbind(seq1,seq2,seq3)
compaired <- list(c("hydra","a_digitifera"),c("Desmophyllum", "a_digitifera"),c("Desmophyllum","hydra"))

#protein heterodimerization activity---
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- data.frame(con_raw1[grepl(con_GO[1,],con_raw1$V4), ])
con_GO1 
GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
# 将序列和对应的名字逐一添加到combined_seqs中
library(seqinr)
combined_data <- c(GO_seq1$seq, GO_seq2$seq, GO_seq3$seq)
names <- c(GO_seq1$ID, GO_seq2$ID, GO_seq3$ID)

write.fasta(combined_data, names=names, file ='coral.fasta',open="a")

combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
ggsave("heterodimerization.pdf",p1, w=6, h=8)
#comparison with transcriptome and proteome-------------
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
extracellular_pI<- seq1[seq1$ID2 %in% unlist(GO1_list_1),]
extracellular_pI$Nr.ID<- extracellular_pI$ID2
extracellular=datahgraw_NR[datahgraw_NR$Nr.ID %in% unlist(GO1_list_1),]
extracellular$geneID=rownames(extracellular)
extracellular_merge <- merge(extracellular_pI,extracellular, by="Nr.ID",all.x =F)
extracellular_merge1 <- extracellular_merge[,c(1,5,12,42)]
normalized_all_d<- as.data.frame(normalized_all)
normalized_all_d$geneID <- rownames(normalized_all_d)
extracellular_heatmap<-merge(normalized_all_d, extracellular_merge1, by="geneID",all.x =F)
extracellular_heatmap$PI1 <-round(extracellular_heatmap$PI1,1)

#RNA CONUT-----
normalized_all$geneID <- rownames(normalized_all)
extracellular_data_merge<- merge(normalized_all,datahgraw,by="geneID")
extracellular_cc <- extracellular_data_merge[grepl("histone H2A",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()

pheatmap(
  extracellular_heatmap[, c(2:5,10:21)],
  scale = 'row',
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = T,
  legend = FALSE,
  labels_row =extracellular_heatmap$PI1,
  color = colorRampPalette(c("blue", "white", "red"))(1000),
)

#extracellular region-----
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- data.frame(con_raw1[grepl(con_GO[2,],con_raw1$V4), ])


GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）

#extracellular all--------------
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- rbind(data.frame(con_raw1[grepl(con_GO[2,],con_raw1$V4), ]),data.frame(ex_raw1[grepl(con_GO[2,],ex_raw1$V4), ]))
GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
p2<- ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
ggsave("extracellular.pdf",p2, w=6, h=8)

##RNA CONUT-----
extracellular_pI<- seq1[seq1$ID2 %in% unlist(GO1_list_1),]
extracellular_pI$Nr.ID<- extracellular_pI$ID2
extracellular=datahgraw_NR[datahgraw_NR$Nr.ID %in% unlist(GO1_list_1),]

extracellular_cc <- extracellular_data_merge[grepl("uncharacterized skeletal organic matrix protein 5",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()
#nucleosome assembly-----
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- data.frame(con_raw1[grepl(con_GO[3,],con_raw1$V4), ])
con_GO1 <- con_GO1
GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")
ALL_list  <- rbind(GO1_list_1,GO1_list_2,GO1_list_3)
GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
compaired <- list(c("hydra","a_digitifera"),c("Desmophyllum", "a_digitifera"),c("Desmophyllum","hydra"))
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）

p3<- ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
ggsave("nuclear.pdf",p3,w=6,h=8)
##RNA CONUT-----
extracellular_pI<- seq1[seq1$ID2 %in% unlist(GO1_list_1),]
extracellular_pI$Nr.ID<- extracellular_pI$ID2
extracellular=datahgraw_NR[datahgraw_NR$Nr.ID %in% unlist(GO1_list_1),]
extracellular
extracellular_cc <- extracellular_data_merge[grepl("histone H2B",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()

# DNA integration"-----
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- data.frame(con_raw1[grepl(con_GO[4,],con_raw1$V4), ])
con_GO1 <- con_GO1
GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
compaired <- list(c("hydra","a_digitifera"),c("Desmophyllum", "a_digitifera"),c("Desmophyllum","hydra"))
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）

ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  theme_classic()
##RNA CONUT-----
extracellular_pI<- seq1[seq1$ID2 %in% unlist(GO1_list_1),]
extracellular_pI<- seq2[seq1$ID2 %in% unlist(GO1_list_2),]
extracellular_pI$Nr.ID<- extracellular_pI$ID2
# max_pi<- extracellular_pI[extracellular_pI$PI1 %in% max(extracellular_pI$PI1),]
# min_pi<-extracellular_pI[extracellular_pI$PI1 %in% min(extracellular_pI$PI1),]

extracellular=datahgraw_NR[datahgraw_NR$Nr.ID %in% unlist(GO1_list_1),]
extracellular
extracellular_cc <- extracellular_data_merge[grepl("uncharacterized protein K02A2.6",extracellular_data_merge$Nr.annotation),]
extracellular_cc$ID2 <- extracellular_cc$Nr.ID
merge(seq1[seq1$ID2 %in% extracellular_cc$Nr.ID,],extracellular_cc, by ="ID2")
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]


# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()



# very long-chain fatty acid biosynthetic process-----
con_GO <- data.frame(genome_con$description[1:10])
con_GO1 <- data.frame(con_raw1[grepl(con_GO[5,],con_raw1$V4), ])
con_GO1 <- con_GO1
GO1_list_1 <- stringr::str_extract_all(con_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(con_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)

anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）

p <- ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(con_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  theme_classic()



#expanded------------
# transepithelial L-ascorbic acid transport"----
ex_GO <- data.frame(genome_ex$description[1:10])
ex_GO1 <- data.frame(ex_raw1[grepl(ex_GO[1,],ex_raw1$V4), ])

GO1_list_1 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
p1<-ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(ex_GO[1,]),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
ggsave("transepithelial.pdf",p1, w=6, h=8)


#RNA CONUT-----

extracellular_cc <- extracellular_data_merge[grepl("solute carrier family 23 member 2",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()


#oxidoreductase activity----
ex_GO <- data.frame(genome_ex$description[1:10])
ex_GO1 <- data.frame(ex_raw1[grepl(ex_GO[2,],ex_raw1$V4), ])

GO1_list_1 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
p1<-ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(ex_GO[2,]),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  ylim(4,16)+
  theme_classic()
ggsave("oxidoreductase.pdf",p1, w=6, h=8)

# Create a bar plot with average and SE----
summary_data <- combined_data %>%
  group_by(group) %>%
  summarise(Average = mean(value),
            SE = sd(value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))

ggplot(summary_data, aes(x = group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = 1, size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()

#RNA CONUT-----

extracellular_cc <- extracellular_data_merge[grepl("dapdiamide synthesis protein DdaC",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()

#indoleacetamide hydrolase activity-----
ex_GO <- data.frame(genome_ex$description[1:10])
ex_GO1 <- data.frame(ex_raw1[grepl(ex_GO[3,],ex_raw1$V4), ])

GO1_list_1 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(ex_GO[3,]),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  geom_hline(aes(yintercept=7))+ 
  ylim(4,16)+
  theme_classic()
ggsave("indoleacetamide.pdf",p2, w=6, h=8)

#RNA CONUT-----
normalized_all$geneID <- rownames(normalized_all)
extracellular_data_merge<- merge(normalized_all,datahgraw,by="geneID")
extracellular_cc <- extracellular_data_merge[grepl("amidase-like",extracellular_data_merge$Nr.annotation),]
extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group
summary_data <- gathered_data %>%
  group_by(Group) %>%
  summarise(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))
compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()


#cell communication"-----
ex_GO <- data.frame(genome_ex$description[1:10])
ex_GO1 <- data.frame(ex_raw1[grepl(ex_GO[4,],ex_raw1$V4), ])

GO1_list_1 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(ex_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  theme_classic()

##RNA CONUT-----
extracellular_pI<- seq1[seq1$ID2 %in% unlist(GO1_list_1),]
extracellular_pI$Nr.ID<- extracellular_pI$ID2
extracellular=datahgraw_NR[datahgraw_NR$Nr.ID %in% unlist(GO1_list_1),]
extracellular
extracellular_cc <- extracellular_data_merge[grepl("ectonucleotide pyrophosphatase/phosphodiesterase family member 5",extracellular_data_merge$Nr.annotation),]

extracellular_cc$row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
row_sums <- rowSums(extracellular_cc[, c(2:5, 14:21)])
# Subset the data frame to keep rows with row sums greater than or equal to 10
filtered_extracellular_ann <- extracellular_cc[extracellular_cc$row_sums %in% max(row_sums), ]#select the highest expressed gene
boxplot_data<- filtered_extracellular_ann[, c(2:5, 14:21)]
# Define the column groupings
con_columns <- c("cont_1_count", "cont_2_count", "cont_3_count", "cont_4_count")
HCHA_columns <- c("CO2_CaOH_1_count", "CO2_CaOH_2_count", "CO2_CaOH_3_count", "CO2_CaOH_4_count")
HCHA1_columns <- c("CO2_NaOH_1_count", "CO2_NaOH_2_count", "CO2_NaOH_3_count", "CO2_NaOH_4_count")
gathered_data <- boxplot_data %>%
  gather(key = "Group", value = "Value", con_columns, HCHA_columns, HCHA1_columns) %>%
  mutate(Group = case_when(
    Group %in% con_columns ~ "con",
    Group %in% HCHA_columns ~ "HCHA",
    Group %in% HCHA1_columns ~ "HCHA1"
  ))
# Create boxplots for each group
anova_result <- aov(Value ~ Group, data = gathered_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
# Calculate the average and SE for each group

compaire_group <- list(c("HCHA","con"),c("HCHA1", "con"))
# Create a bar plot with average and SE
ggplot(summary_data, aes(x = Group, y = Average)) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(position = position_dodge(width = 0.8), shape = "-", size = 9) +  # Short lines for averages
  geom_signif(comparisons = compaire_group,step_increase = 0.2,map_signif_level = T)+
  theme_classic()




#protein folding
ex_GO <- data.frame(genome_ex$description[1:10])
ex_GO1 <- data.frame(ex_raw1[grepl(ex_GO[5,],ex_raw1$V4), ])

GO1_list_1 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=a_digitifera\\|)([^|;]+)")
GO1_list_2 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Hydra_vulgaris_fresh_water_hydra\\|)([^|;]+)")
GO1_list_3 <- stringr::str_extract_all(ex_GO1$V5,  "(?<=Desmophyllum_pertusum_deepsea\\|)([^|;]+)")

GO_seq1<-  seq1[seq1$ID2 %in% unlist(GO1_list_1),]
GO_seq2<- seq2[seq2$ID2 %in% unlist(GO1_list_2),]
GO_seq3<- seq3[seq3$ID2 %in% unlist(GO1_list_3),]
combined_data <- rbind(
  data.frame(value = GO_seq1$PI1, group = "a_digitifera", ID=GO_seq1$ID2),
  data.frame(value = GO_seq2$PI1, group = "hydra", ID=GO_seq2$ID2),
  data.frame(value = GO_seq3$PI1, group = "Desmophyllum", ID=GO_seq3$ID2)
)
anova_result <- aov(value ~ group, data = combined_data)# 运行单因素方差分析（ANOVA）
tukey_result <- TukeyHSD(anova_result)# 进行事后多重比较校正（Tukey's HSD）
ggplot(combined_data, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point() +
  labs(title = paste0(ex_GO1$V4),
       x = "species",
       y = "pI Value")+
  # geom_text(aes(label = ID), hjust = - 0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.2,map_signif_level = T)+
  theme_classic()

#freshwater expanded-----

genome_ex_ALL <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/genome_GO.csv', header = T) 
genome_ex_ALL$description <- trimws(genome_ex_ALL$description)
genome_ex_ALL_e <- genome_ex_ALL[genome_ex_ALL$habitat=="freshwater" & genome_ex_ALL$contracted.expanded=="expanded",]
list_species <- unique(genome_ex_ALL_e$species)
list_1 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[1],]$description
list_2 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[2],]$description
list_3 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[3],]$description
list_4 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[4],]$description
list_5 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[5],]$description
list_6 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[6],]$description


# 将六个列表存储到一个列表中
all_lists <- list(sponge=list_3, coral=list_4, snail=list_5)
inter_lists <- get.venn.partitions(all_lists)
combined_list<- inter_lists[inter_lists$sponge=="TRUE"&inter_lists$coral=="TRUE"&inter_lists$snail=="TRUE",]$..values..
combined_list1 <- data.frame(des = combined_list)
venn_plot <-ggvenn(all_lists,
                   show_percentage =F)
combined_list1$X1
venn_plot <- venn_plot +
  annotate(geom = "text", x = 0.5, y = 0.5, label = combined_list1$X1, color = "red", size = 4)
ggsave("fresh_expanded.pdf",venn_plot )
#freshwater contracted-----
genome_ex_ALL_e <- genome_ex_ALL[genome_ex_ALL$habitat=="freshwater" & genome_ex_ALL$contracted.expanded=="contracted",]

list_species <- unique(genome_ex_ALL_e$species)
list_1 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[1],]$description
list_2 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[2],]$description
list_3 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[3],]$description
list_4 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[4],]$description
list_5 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[5],]$description
list_6 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[6],]$description


# 将六个列表存储到一个列表中
all_lists <- list(sponge=list_3, coral=list_4, snail=list_5)
inter_lists <- get.venn.partitions(all_lists)
combined_list<- inter_lists[inter_lists$sponge=="TRUE"&inter_lists$coral=="TRUE"&inter_lists$snail=="TRUE",]$..values..
combined_list1 <- data.frame(des = combined_list)
venn_plot <-ggvenn(all_lists,
                   show_percentage =F)
combined_list1$X1
venn_plot <- venn_plot +
  annotate(geom = "text", x = 0.5, y = 0.5, label = combined_list1$X1, color = "red", size = 4)
ggsave("fresh_contracted.pdf",venn_plot )

#deepsea expanded-----
genome_ex_ALL_e <- genome_ex_ALL[genome_ex_ALL$habitat=="deepsea" & genome_ex_ALL$contracted.expanded=="expanded",]

list_species <- unique(genome_ex_ALL_e$species)
list_1 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[1],]$description
list_2 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[2],]$description
list_3 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[3],]$description
list_4 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[4],]$description

# 将六个列表存储到一个列表中
all_lists <- list(sponge=list_2, coral=list_3, snail=list_4)
inter_lists <- get.venn.partitions(all_lists)
combined_list<- inter_lists[inter_lists$sponge=="TRUE"&inter_lists$coral=="TRUE"&inter_lists$snail=="TRUE",]$..values..
combined_list1 <- data.frame(des = combined_list)
venn_plot <-ggvenn(all_lists,
                   show_percentage =F)

combined_list1$X1
venn_plot <- venn_plot +
  annotate(geom = "text", x = 0.5, y = 0.5, label = combined_list1$X1, color = "red", size = 4)
ggsave("deepsea_expanded.pdf",venn_plot )


#deepsea contracted-----
genome_ex_ALL_e <- genome_ex_ALL[genome_ex_ALL$habitat=="deepsea" & genome_ex_ALL$contracted.expanded=="contracted",]

list_species <- unique(genome_ex_ALL_e$species)
list_2 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[2],]$description
list_3 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[3],]$description
list_4 <- genome_ex_ALL_e[genome_ex_ALL_e$species==list_species[4],]$description


# 将六个列表存储到一个列表中
all_lists <- list(sponge=list_2, coral=list_3, snail=list_4)
inter_lists <- get.venn.partitions(all_lists)
combined_list<- inter_lists[inter_lists$sponge=="TRUE"&inter_lists$coral=="TRUE"&inter_lists$snail=="TRUE",]$..values..
combined_list1 <- data.frame(des = combined_list)
venn_plot <-ggvenn(all_lists,
                   show_percentage =F)
combined_list1$X1
venn_plot <- venn_plot +
  annotate(geom = "text", x = 0.5, y = 0.5, label = combined_list1$X1, color = "red", size = 4)
ggsave("deepsea_contracted.pdf",venn_plot )


