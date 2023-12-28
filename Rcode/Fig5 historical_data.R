library(ggplot2)
library("Biostrings")
library(Peptides) 
library(rentrez)
library(pheatmap)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
# 定义处理的条件
conditions <- c("CaOH2_up", "CaOH2_down", "CO2_up", "CO2_down", "CO2CaOH2_up", "CO2CaOH2_down","CO2NaOH_up", "CO2NaOH_down")

# 循环遍历每个条件
for (condition in conditions) {
  # 获取相应条件的Accession号
  condition_name <- get(condition)
  accession_numbers <- datahgraw[datahgraw$geneID %in% condition_name$name, ]$Nr.ID[datahgraw[datahgraw$geneID %in% condition_name$name, ]$Nr.ID != "--"]
  
  # 初始化字符向量
  protein_sequences <- character(0)
  
  # 循环遪Accession号并下载蛋白序列
  for (accession in accession_numbers) {
    # 使用efetch函数下载蛋白序列
    protein_seq <- entrez_fetch(db = "protein", id = accession, rettype = "fasta", retmode = "text")
    # 打印蛋白序列
    cat(protein_seq, "\n")
    # 将蛋白序列附加到字符向量中
    protein_sequences <- c(protein_sequences, protein_seq)
  }
  
  # 保存蛋白序列到相应的文件
  write(paste(protein_sequences, collapse = "\n"), file = paste0(condition, "_protein_sequences.fasta"))
}

#CO2 fasta----
fastaFile_CO2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi -----
seq_gene_CO2_up=data.frame(pi=pI(fastaFile_CO2_up))
seq_gene_CO2_up$names=  substr(names(fastaFile_CO2_up),1,14)

seq_gene_CO2_down=data.frame(pi=pI(fastaFile_CO2_down))
 
ggplot() +
  geom_density(data = seq_gene_CO2_up, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = seq_gene_CO2_down, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))





#CaOH2----
fastaFile_CaOH2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CaOH2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CaOH2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CaOH2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
seq_gene_CaOH2_up=data.frame(pi=pI(fastaFile_CaOH2_up))
seq_gene_CaOH2_down=data.frame(pi=pI(fastaFile_CaOH2_down))
ggplot() +
  geom_density(data = seq_gene_CaOH2_up, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = seq_gene_CaOH2_down, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))
#CO2CaOH2----
fastaFile_CO2CaOH2_up <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2CaOH2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
fastaFile_CO2CaOH2_down <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2CaOH2_down_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
seq_gene_CO2CaOH2_up=data.frame(pi=pI(fastaFile_CO2CaOH2_up))
seq_gene_CO2CaOH2_down=data.frame(pi=pI(fastaFile_CO2CaOH2_down))
ggplot() +
  geom_density(data = seq_gene_CO2CaOH2_up, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = seq_gene_CO2CaOH2_down, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))


#data from: sun-----
#"Sun, Y., L. Jiang, S. Gong, G. Diaz-Pulido, X. Yuan, H. Tong, L. Huang, G. Zhou, Y. Zhang & H. Huang (2022) Changes in physiological performance and protein expression in the larvae of the coral Pocillopora damicornis and their symbionts in response to elevated temperature and acidification. Sci Total Environ, 807, 151251."
raw_data <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/sun_data.csv', header = T) 
head(raw_data)
up_gene <- raw_data[raw_data$Fold.change.OA.Control>1.2,]
down_gene <- raw_data[raw_data$Fold.change.OA.Control<0.8,]
# 定义Accession号码
accession_numbers <- raw_data$Accession
protein_sequences <- character(0)
# 循环遍历Accession号并下载蛋白序列
for (accession in accession_numbers) {
  # 使用efetch函数下载蛋白序列
  protein_seq <- entrez_fetch(db="protein", id=accession, rettype="fasta", retmode="text")
  # 打印蛋白序列
  cat(protein_seq, "\n")
  # 将蛋白序列附加到字符向量中
  protein_sequences <- c(protein_sequences, protein_seq)
}

# 将所有蛋白序列保存到一个文件
write(paste(protein_sequences, collapse = "\n"), file = "all_protein_sequences.fasta")
fastaFile_CO2 <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/CO2_up_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
seq_gene=data.frame(pi=pI(fastaFile1))
seq_gene$names=substr(fastaFile1@ranges@NAMES, 1, 14)
up_seq_gene <- seq_gene[seq_gene$names %in% up_gene$Accession,]
down_seq_gene <- seq_gene[seq_gene$names %in% down_gene$Accession,]
ggplot() +
  geom_density(data = up_seq_gene, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = down_seq_gene, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))
#seqs <- readDNAStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/Galaxea_Unigene.fa', format = "fasta", use.names = TRUE)
fastaFile1 <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/all_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
seq_gene=data.frame(pi=pI(fastaFile1))
seq_gene$names=substr(fastaFile1@ranges@NAMES, 1, 14)
up_seq_gene <- seq_gene[seq_gene$names %in% up_gene$Accession,]
down_seq_gene <- seq_gene[seq_gene$names %in% down_gene$Accession,]
ggplot() +
  geom_density(data = up_seq_gene, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = down_seq_gene, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))
#data from: lin-----
# data from "Lin, Z., L. Wang, M. Chen, X. Zheng & J. Chen (2022) Proteome and microbiota analyses characterizing dynamic coral-algae-microbe tripartite interactions under simulated rapid ocean acidification. Sci Total Environ, 810, 152266."

raw_data2 <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/lin_data.csv', header = T) 
head(raw_data2$Protein.description)
pattern <- "gb\\|([^|]+)\\|"
matches <- regmatches(raw_data2$Protein.description, gregexpr(pattern, raw_data2$Protein.description))
result_ID <- gsub(pattern, "\\1", matches)
raw_data2$ID  <-result_ID
pattern <- "^E.*"
filtered_result_ID<- result[grep(pattern, result_ID)]


accession_numbers <- filtered_result_ID
protein_sequences <- character(0)

# 循环遍历Accession号并下载蛋白序列
for (accession in accession_numbers) {
  # 使用efetch函数下载蛋白序列
  protein_seq <- entrez_fetch(db="protein", id=accession, rettype="fasta", retmode="text")
  # 打印蛋白序列
  cat(protein_seq, "\n")
  # 将蛋白序列附加到字符向量中
  protein_sequences <- c(protein_sequences, protein_seq)
}

# 将所有蛋白序列保存到一个文件
write(paste(protein_sequences, collapse = "\n"), file = "lin_protein_sequences.fasta")

#seqs <- readDNAStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/Galaxea_Unigene.fa', format = "fasta", use.names = TRUE)
fastaFile2 <- readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/lin_protein_sequences.fasta', format = "fasta", use.names = TRUE)
#pi of genome, mRNA and protein
up_gene2 <- raw_data2[raw_data2$fold1>1.2,]
down_gene2 <- raw_data2[raw_data2$fold1<0.8,]
seq_gene=data.frame(pi=pI(fastaFile2))
seq_gene$names=substr(fastaFile2@ranges@NAMES, 1, 10)
up_seq_gene <- seq_gene[seq_gene$names %in% up_gene2$ID,]
down_seq_gene <- seq_gene[seq_gene$names %in% down_gene2$ID,]
P1<- ggplot() +
  geom_density(data = up_seq_gene, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = down_seq_gene, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+scale_y_continuous(limits = c(0, 0.3))+
  labs(subtitle = "Data from Lin et al. 2022, pH = 7.85 vs pH = 8.15")


#pi of genome, mRNA and protein
up_gene2 <- raw_data2[raw_data2$fold2>1.2,]
down_gene2 <- raw_data2[raw_data2$fold2<0.8,]
seq_gene=data.frame(pi=pI(fastaFile2))
seq_gene$names=substr(fastaFile2@ranges@NAMES, 1, 10)
up_seq_gene <- seq_gene[seq_gene$names %in% up_gene2$ID,]
down_seq_gene <- seq_gene[seq_gene$names %in% down_gene2$ID,]
P2<- ggplot() +
  geom_density(data = up_seq_gene, aes(x = pi, colour = "upregulated"), alpha = 2, size = 1.5) +
  geom_density(data = down_seq_gene, aes(x = pi, colour = "downregulated"), alpha = 2, size = 1.5)+
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  theme_classic()+
  scale_x_continuous(limits = c(4, 12))+
  scale_y_continuous(limits = c(0, 0.3))+
  labs(subtitle = "Data from Lin et al. 2022, pH = 7.85 vs pH = 8.45")
combine_plot<- P1|P2
ggsave("lindata.pdf",combine_plot,w=8,h=6)
