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
# mollusca organism------
a <- list.files("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca")  
a
fastaFile1<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca/Biomphalaria glabrata seawater mollusca sea snail.fasta', format = "fasta", use.names = TRUE)
fastaFile2<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca/Gigantopelta aegis deep sea vent snail.fasta', format = "fasta", use.names = TRUE)
fastaFile3<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca/Haliotis rufescens seawater snail.fasta', format = "fasta", use.names = TRUE)
fastaFile4<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca/Lottia gigantea seawater snail.fasta', format = "fasta", use.names = TRUE)
fastaFile5<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/mollusca/Pomacea canaliculata freshwater snail.fasta', format = "fasta", use.names = TRUE)



name1 <- names(fastaFile1) %>% as.data.frame()
name2 <- names(fastaFile2) %>% as.data.frame()
name3 <- names(fastaFile3) %>% as.data.frame()
name4 <- names(fastaFile4) %>% as.data.frame()
name5 <- names(fastaFile5) %>% as.data.frame()


data1 <- as.data.frame(fastaFile1)
data2 <- as.data.frame(fastaFile2)
data3 <- as.data.frame(fastaFile3)
data4 <- as.data.frame(fastaFile4)
data5 <- as.data.frame(fastaFile5)


data21 <- data.frame(name1,seq=data1$x)
data22 <- data.frame(name2,seq=data2$x)
data23 <- data.frame(name3,seq=data3$x)
data24 <- data.frame(name4,seq=data4$x)
data25 <- data.frame(name5,seq=data5$x)



seq1=data.frame(data21) %>%
  mutate("PI" = pI(seq = data21$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data21$seq))
seq2=data.frame(data22) %>%
  mutate("PI" = pI(seq = data22$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data22$seq))
seq3=data.frame(data23) %>%
  mutate("PI" = pI(seq = data23$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data23$seq))
seq4=data.frame(data24) %>%
  mutate("PI" = pI(seq = data24$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data24$seq))
seq5=data.frame(data25) %>%
  mutate("PI" = pI(seq = data25$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data25$seq))


library("plyr") 
list1 <- list()
list1[[1]] <- data.frame(t(seq1$PI))
list1[[2]] <- data.frame(t(seq2$PI))
list1[[3]] <- data.frame(t(seq3$PI))
list1[[4]] <- data.frame(t(seq4$PI))
list1[[5]] <- data.frame(t(seq5$PI))

u <- rbind.fill(list1)
data_all=data.frame(t(u))
head(data_all)
data_all_m=colMeans(na.omit(data_all))
boxplot(data_all)


a
pdf("density_snail.pdf",width=15,height=6)

ggplot(data_all)+ 
  geom_density(aes(x=X1,colour = "Biomphalaria glabrata seawater mollusca sea snail"),alpha = 2,size=1.5) +
  geom_density(aes(x=X2,colour = "Gigantopelta aegis deep sea vent snail"),alpha = 2, linetype = 2,size=1) +
  geom_density(aes(x=X3,colour = "Haliotis rufescens seawater snail"),alpha = 2,size=1.5) +
  geom_density(aes(x=X4,colour = "Lottia gigantea seawater snail"),alpha = 2,size=1.5) +
  geom_density(aes(x=X5,colour = "Pomacea canaliculata freshwater snail"),alpha = 2, linetype = 2,size=1) +
  
  #scale_color_brewer(palette = "Paired", direction = -1)+ 
  scale_x_continuous(breaks=seq(0,15,1))+
  scale_y_continuous(breaks=seq(0,0.38,0.05),limits = c(0,0.38))+
  theme_classic()

dev.off() 




