library(tidyverse)
library(Biostrings)  # 用于读取fasta格式
library(Peptides)  
#Osorio, D., Rondon-Villarreal, P. & Torres, R. Peptides: A package for data mining of antimicrobial peptides. The R Journal. 7(1), 4–14 (2015).
library(seqinr)
library(stringi)
library(seqinr)
library(ggplot2)
library(ggsignif)
library(vioplot)
library(WGCNA)
library("Biostrings")
library(pheatmap) 
library(RColorBrewer)
library(ggcorrplot)
library(grid)
library(ape)
library(ggseqlogo)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence')
# sponge organism------
a <- list.files("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/sponge")  
a
fastaFile1<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/sponge/Amphimedon queenslandica intertidal sea.fasta', format = "fasta", use.names = TRUE)
fastaFile2<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/sponge/Oopsacas minuta deep sea glass sponge.fasta', format = "fasta", use.names = TRUE)

name1 <- names(fastaFile1) %>% as.data.frame()
name2 <- names(fastaFile2) %>% as.data.frame()


data1 <- as.data.frame(fastaFile1)
data2 <- as.data.frame(fastaFile2)

data21 <- data.frame(name1,seq=data1$x)
data22 <- data.frame(name2,seq=data2$x)

seq1=data.frame(data21) %>%
  mutate("PI" = pI(seq = data21$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data21$seq))
seq2=data.frame(data22) %>%
  mutate("PI" = pI(seq = data22$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data22$seq))

library("plyr") 
list1 <- list()
list1[[1]] <- data.frame(t(seq1$PI))
list1[[2]] <- data.frame(t(seq2$PI))

u <- rbind.fill(list1)
data_all=data.frame(t(u))
head(data_all)
data_all_m=colMeans(na.omit(data_all))
boxplot(data_all)

pdf("density_sponge.pdf",width=15,height=6)
a
ggplot(data_all)+ 
  geom_density(aes(x=X1,colour = "Amphimedon queenslandica intertidal sea.fasta"),alpha = 2,size=1.5) +
  geom_density(aes(x=X2,colour = "Oopsacas minuta deep sea glass sponge.fasta"),alpha = 2,size=1.5) +
  scale_x_continuous(breaks=seq(0,15,1))+
  scale_y_continuous(breaks=seq(0,0.38,0.05),limits = c(0,0.38))+
  theme_classic()
dev.off() 




