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
# cnidarian organism------
a <- list.files("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian")  
a
fastaFile1<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Acropora millepora sea coral.fasta', format = "fasta", use.names = TRUE)
fastaFile2<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Pocillopora damicornis seawater coral.fasta', format = "fasta", use.names = TRUE)
fastaFile3<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Stylophora pistillata seawater coral.fasta', format = "fasta", use.names = TRUE)
fastaFile4<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Dendronephthya gigantea soft coral 10-20m.fasta', format = "fasta", use.names = TRUE)
fastaFile5<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Nematostella vectensis sea anemone.fasta', format = "fasta", use.names = TRUE)
fastaFile6<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Hydra vulgaris fresh water hydra.fasta', format = "fasta", use.names = TRUE)
fastaFile7<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/Hydra viridissima freshwater hydra.fasta', format = "fasta", use.names = TRUE)


name1 <- names(fastaFile1) %>% as.data.frame()
name2 <- names(fastaFile2) %>% as.data.frame()
name3 <- names(fastaFile3) %>% as.data.frame()
name4 <- names(fastaFile4) %>% as.data.frame()
name5 <- names(fastaFile5) %>% as.data.frame()
name6 <- names(fastaFile6) %>% as.data.frame()
name7 <- names(fastaFile7) %>% as.data.frame()

data1 <- as.data.frame(fastaFile1)
data2 <- as.data.frame(fastaFile2)
data3 <- as.data.frame(fastaFile3)
data4 <- as.data.frame(fastaFile4)
data5 <- as.data.frame(fastaFile5)
data6 <- as.data.frame(fastaFile6)
data7 <- as.data.frame(fastaFile7)

data21 <- data.frame(name1,seq=data1$x)
data22 <- data.frame(name2,seq=data2$x)
data23 <- data.frame(name3,seq=data3$x)
data24 <- data.frame(name4,seq=data4$x)
data25 <- data.frame(name5,seq=data5$x)
data26 <- data.frame(name6,seq=data6$x)
data27 <- data.frame(name7,seq=data7$x)


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
seq6=data.frame(data26) %>%
  mutate("PI" = pI(seq = data26$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data26$seq))
seq7=data.frame(data27) %>%
  mutate("PI" = pI(seq = data27$seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = data27$seq))


library("plyr") 
list1 <- list()
list1[[1]] <- data.frame(t(seq1$PI))
list1[[2]] <- data.frame(t(seq2$PI))
list1[[3]] <- data.frame(t(seq3$PI))
list1[[4]] <- data.frame(t(seq4$PI))
list1[[5]] <- data.frame(t(seq5$PI))
list1[[6]] <- data.frame(t(seq6$PI))
list1[[7]] <- data.frame(t(seq7$PI))
u <- rbind.fill(list1)
data_all=data.frame(t(u))
head(data_all)
data_all_m=colMeans(na.omit(data_all))
boxplot(data_all)
a
pdf("density_coral_hydra.pdf",width=15,height=6)

ggplot(data_all)+ 
  geom_density(aes(x=X1,colour = "Acropora millepora,stony coral"),alpha = 2,size=1.5) +
  geom_density(aes(x=X2,colour = "Pocillopora damicornis,stony coral"),alpha = 2,size=1.5) +
  geom_density(aes(x=X3,colour = "Stylophora pistillata,stony coral"),alpha = 2,size=1.5) +
  geom_density(aes(x=X4,colour = "Dendronephthya gigantea, soft coral"),alpha = 2,size=1.5) +
  geom_density(aes(x=X5,colour = "Nematostella vectensis, sea anemone"),alpha = 2,size=1.5) +
  geom_density(aes(x=X7,colour = "Hydra viridissima freshwater hydra"),alpha = 2,linetype = 2,size=1) +
  geom_density(aes(x=X6,colour = "Hydra vulgaris, fresh water hydra"),alpha = 2, linetype = 2,size=1) +
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  scale_x_continuous(breaks=seq(0,15,1))+
  scale_y_continuous(breaks=seq(0,0.38,0.05),limits = c(0,0.38))+
  theme_classic()
dev.off() 








# allnames <-c("Acropora millepora",
#              'Pocillopora damicornis',
#              'Stylophora pistillata',
#              'Nematostella vectensis',
#              'Crassostrea gigas',
#              'Amphiprion percula',
#              'Amphimedon queenslandica',
#              'Ectocarpus siliculosus',
#              'Helicobacter pylori',
#              'Symbiodinium necroappetens')
# #plot charge----
# seq1_<-seq1[which(seq1$charge<50&seq1$charge>-50),4]
# seq2_<-seq2[which(seq2$charge<50&seq2$charge>-50),4]
# seq3_<-seq3[which(seq3$charge<50&seq3$charge>-50),4]
# seq4_<-seq4[which(seq4$charge<50&seq4$charge>-50),4]
# seq5_<-seq5[which(seq5$charge<50&seq5$charge>-50),4]
# seq6_<-seq6[which(seq6$charge<50&seq6$charge>-50),4]
# seq7_<-seq7[which(seq7$charge<50&seq7$charge>-50),4]
# seq8_<-seq8[which(seq8$charge<50&seq8$charge>-50),4]
# seq9_<-seq9[which(seq9$charge<50&seq9$charge>-50),4]
# seq10_<-seq10[which(seq10$charge<50&seq10$charge>-50),4]
# 
# 
# list2 <- list()
# list2[[1]] <- data.frame(t(seq1_))
# list2[[2]] <- data.frame(t(seq2_))
# list2[[3]] <- data.frame(t(seq3_))
# list2[[4]] <- data.frame(t(seq4_))
# list2[[5]] <- data.frame(t(seq5_))
# list2[[6]] <- data.frame(t(seq6_))
# list2[[7]] <- data.frame(t(seq7_))
# list2[[8]] <- data.frame(t(seq8_))
# list2[[9]] <- data.frame(t(seq9_))
# list2[[10]] <- data.frame(t(seq10_))
# 
# u2 <- rbind.fill(list2)
# data_all2=data.frame(t(u2))
# head(data_all2)
# allnames <-c("Acropora millepora",
#              'Pocillopora damicornis',
#              'Stylophora pistillata',
#              'Nematostella vectensis',
#              'Crassostrea gigas',
#              'Amphiprion percula',
#              'Amphimedon queenslandica',
#              'Ectocarpus siliculosus',
#              'Helicobacter pylori',
#              'Symbiodinium necroappetens')
# 
# ggplot(data_all2)+ 
#   geom_density(aes(x=X1,colour = "Acropora millepora"),alpha = 3,size=2) +
#   geom_density(aes(x=X2,colour = "Pocillopora damicornis"),alpha = 2,size=2) +
#   geom_density(aes(x=X3,colour = "Stylophora pistillata"),alpha = 2,size=2) +
#   geom_density(aes(x=X4,colour = "Nematostella vectensis"),alpha = 2) +
#   geom_density(aes(x=X5,colour = "Crassostrea gigas"),alpha = 2) +
#   geom_density(aes(x=X6,colour = "Amphiprion percula"),alpha = 2) +
#   #geom_density(aes(x=X7,colour = "Amphimedon queenslandica"),alpha = 2) +
#   geom_density(aes(x=X8,colour = "Ectocarpus siliculosus"),alpha = 2,linetype = 2,size=1) +
#   #geom_density(aes(x=X9,colour = "Helicobacter pylori"),alpha = 2,linetype = 2,size=1) +
#   geom_density(aes(x=X10,colour = "Symbiodinium necroappetens"),alpha = 2, linetype = 2,size=1) +
#   theme_classic()


# fastaFile<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/Ectocarpus siliculosus protein.fasta', format = "fasta", use.names = TRUE)
# name <- names(fastaFile) %>% as.data.frame()
# data1 <- as.data.frame(fastaFile)
# data2 <- data.frame(name=name$.,seq=data1$x)
# data_rownames=data2[grep(pattern="tRNA ligase",data2$name),]
# seqPI=data.frame(data_rownames) %>%
#   mutate("length" = Peptides::lengthpep(seq = data_rownames$seq)) %>%  # lengthpep() 计算长度
#   mutate("mw" = mw(seq = data_rownames$seq)) %>%                 # mw() 计算分子量
#   mutate("hydrophobicity" = hydrophobicity(seq =  data_rownames$seq)) %>%     # hydrophobicity() 计算疏水性
#   mutate("PI" = pI(seq = data_rownames$seq)) %>%                 # pI() 计算等电点
#   mutate("charge" = charge(seq = data_rownames$seq))%>% 
#   mutate("instaIndex" = instaIndex(seq =  data_rownames$seq))
# ggplot(seqPI, aes(x=PI))+ 
#   geom_density(alpha = 2) +
#   theme_classic() 
