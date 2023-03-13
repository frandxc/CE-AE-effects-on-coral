#BiocManager::install("PerformanceAnalytics")
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
library(ggplot2)
library(ggsignif)
library(vioplot)
library(WGCNA)
library("Biostrings")
library(pheatmap) 
library(RColorBrewer)
library(ggcorrplot)
library(Peptides) 
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence')

#all_fasta<- read.fasta('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/Galaxea_Unigene.fa')


load('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/gene and protein correlation/corhghp1.RData')
load('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/heatmap/heatmap.RData')
#seqs <- readDNAStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/raw data/Galaxea_Unigene.fa', format = "fasta", use.names = TRUE)
fastaFile1<- Biostrings::readBStringSet('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R/sequence/cnidarian/a_digitifera.fasta', format = "fasta", use.names = TRUE)
#seqs_translated <- translate(seqs)# translate mRNA into protein

#选取
# corhp_Alkalinity1=rbind(corhp_Alkalinity[1:10,],corhp_Alkalinity[996:1005,])
seq_gene=data.frame(pi=pI(fastaFile1))
row.names(seq_gene)=substr(row.names(seq_gene), 1, 15)

aa=as.character(row.names(datahp_log_cen) )
p_Nr=datahgraw1[aa,]$Nr.ID
mRNA_Nr=datahgraw1$Nr.ID
seq_p=data.frame(pi=c(seq_gene[p_Nr,]))
seq_mRNA=data.frame(pi=c(seq_gene[mRNA_Nr,]))
# letterFrequency(DNA_seq[1],"A")
# paste(DNA_seq[[1]], collapse="") #extract sequence
# nchar(DNA_seq)
library("plyr") 
#ggplot(seq_mRNA)+ 
list1 <- list()
list1[[1]] <- data.frame(t(seq_gene$pi))
list1[[2]] <- data.frame(t(seq_p$pi))
list1[[3]] <- data.frame(t(seq_mRNA$pi))

u <- rbind.fill(list1)
data_all=data.frame(t(u))


pdf("mRNA_protein.pdf",width=10,height=6)
ggplot(data_all)+ 
  geom_density(aes(x=X1,colour = "genome"),alpha = 2,size=1.5) +
  geom_density(aes(x=X2,colour = "Protein"),alpha = 2,size=1.5) +
  geom_density(aes(x=X3,colour = "mRNA"),alpha = 2,size=1.5) +
  scale_color_brewer(palette = "Paired", direction = -1)+ 
  scale_x_continuous(breaks=seq(0,15,2))+
  scale_y_continuous(breaks=seq(0,0.23,0.05),limits = c(0,0.23))+
  theme_classic()
dev.off() 


# protein cor with traits

compaired <- list(c("Control", "Ca(OH)2"), 
                  c("Control","CO2"), 
                  c("Control","CO2+Ca(OH)2"),
                  c("Control","CO2+NaOH"))

row.names(seq_p) <-row.names(datahp_log_cen) 

#protein_data_s=t(scale(t(datahp_log_cen),center=T,scale=F))

pi_protein<-merge(seq_p,datahp_log_cen,by='row.names',all=F)

pi_protein10<-subset(pi_protein,pi>7.5)
pi_protein10_m<-data.frame(mean=colMeans(pi_protein10[,3:22]))
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
pi_protein10_m2 <- data.frame(conditionshg,abundance=pi_protein10_m)
pi_protein10_m2$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pdf("pi7.5.pdf",width=6,height=6) 
ggplot(pi_protein10_m2,aes(conditionshg,mean,color=conditionshg),y="Correlation R2", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
#7-7.5----
pi_protein10<-subset(pi_protein,pi>7&pi<7.5)
pi_protein10_m<-data.frame(mean=colMeans(pi_protein10[,3:22]))
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
pi_protein10_m2 <- data.frame(conditionshg,abundance=pi_protein10_m)
pi_protein10_m2$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pdf("pi7-7.5.pdf",width=6,height=6) 
ggplot(pi_protein10_m2,aes(conditionshg,mean,color=conditionshg),y="Correlation R2", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#<7----
pi_protein10<-subset(pi_protein,pi<7)
pi_protein10_m<-data.frame(mean=colMeans(pi_protein10[,3:22]))
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
pi_protein10_m2 <- data.frame(conditionshg,abundance=pi_protein10_m)
pi_protein10_m2$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pdf("pi7.pdf",width=6,height=6) 
ggplot(pi_protein10_m2,aes(conditionshg,mean,color=conditionshg),y="Correlation R2", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#all----
pi_protein10<-pi_protein
pi_protein10_m<-data.frame(mean=colMeans(pi_protein10[,3:22]))
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
pi_protein10_m2 <- data.frame(conditionshg,abundance=pi_protein10_m)
pi_protein10_m2$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
pdf("all.pdf",width=6,height=6) 
ggplot(pi_protein10_m2,aes(conditionshg,mean,color=conditionshg),y="Correlation R2", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

