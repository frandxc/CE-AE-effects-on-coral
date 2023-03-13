library(tidyverse)
library(Biostrings)  # 用于读取fasta格式
library(Peptides)  
#Osorio, D., Rondon-Villarreal, P. & Torres, R. Peptides: A package for data mining of antimicrobial peptides. The R Journal. 7(1), 4–14 (2015).
library(ggplot2)
library(ggsignif)
library(vioplot)
library(WGCNA)
library("Biostrings")
library(pheatmap) 
library(RColorBrewer)
library(ggcorrplot)
library(grid)
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig3.RData")
AA_seqs <-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/peptide.csv', header = T) 
traits <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/traits.csv', header = T, row.names=1) 
protein_seqs <-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/protein_LFQ.csv', header = T) 
rownames(AA_seqs)<-AA_seqs $Leading.razor.protein
fa=data.frame(name=AA_seqs$Leading.razor.protein,seq=AA_seqs$Sequence)
seqPI=data.frame(fa) %>%
  mutate("length" = Peptides::lengthpep(seq= seq)) %>%  # lengthpep() 计算长度
  mutate("mw" = mw(seq= seq)) %>%                 # mw() 计算分子量
  mutate("hydrophobicity" = hydrophobicity(seq= seq)) %>%     # hydrophobicity() 计算疏水性
  mutate("PI" = pI(seq= seq)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = seq))%>% 
  mutate("instaIndex" = instaIndex(seq= seq))

rownames(seqPI)<-AA_seqs$Leading.razor.protein
normalized_hp_d<-datahp_log_cen
seqpi_abundance <- merge(seqPI,normalized_hp_d,by='row.names',all=F)
seqpi_traits<- data.frame(ID=seqpi_abundance$Row.names,seqpi_abundance[10:29] )


DNA_seq=AA_seqs$Sequence
# letterFrequency(DNA_seq[1],"A")
paste(DNA_seq[[1]], collapse="") #extract sequence
nchar(DNA_seq)
#Unigene0168205
library(stringi)
AG=stri_count_regex(DNA_seq, "G")/nchar(DNA_seq)*100
AA=stri_count_regex(DNA_seq, "A")/nchar(DNA_seq)*100
AV=stri_count_regex(DNA_seq, "V")/nchar(DNA_seq)*100
AL=stri_count_regex(DNA_seq, "L")/nchar(DNA_seq)*100
AI=stri_count_regex(DNA_seq, "I")/nchar(DNA_seq)*100
AF=stri_count_regex(DNA_seq, "F")/nchar(DNA_seq)*100
AW=stri_count_regex(DNA_seq, "W")/nchar(DNA_seq)*100
AY=stri_count_regex(DNA_seq, "Y")/nchar(DNA_seq)*100
AD=stri_count_regex(DNA_seq, "D")/nchar(DNA_seq)*100
AN=stri_count_regex(DNA_seq, "N")/nchar(DNA_seq)*100
AE=stri_count_regex(DNA_seq, "E")/nchar(DNA_seq)*100
AK=stri_count_regex(DNA_seq, "K")/nchar(DNA_seq)*100
AQ=stri_count_regex(DNA_seq, "Q")/nchar(DNA_seq)*100
AM=stri_count_regex(DNA_seq, "M")/nchar(DNA_seq)*100
AS=stri_count_regex(DNA_seq, "S")/nchar(DNA_seq)*100
AT=stri_count_regex(DNA_seq, "T")/nchar(DNA_seq)*100
AC=stri_count_regex(DNA_seq, "C")/nchar(DNA_seq)*100
AP=stri_count_regex(DNA_seq, "P")/nchar(DNA_seq)*100
AH=stri_count_regex(DNA_seq, "H")/nchar(DNA_seq)*100
AR=stri_count_regex(DNA_seq, "R")/nchar(DNA_seq)*100
ami= c("D","E","K","H","R","G","A","V","L","I","M","F","W","P","S","T","C","Y","N","Q")

value=data.frame(AD,AE,AK,AH,AR,AG,AA,AV,AL,AI,AM,AF,AW,AP,AS,AT,AC,AY,AN,AQ)
row.names(value)=AA_seqs$Leading.razor.protein
aa_PI=merge(value,seqPI,by=0,all =F)
rownames(aa_PI)=aa_PI$Row.names
conditionshg <- factor(c(rep("Control",4),rep( "Ca(OH)2",4),rep("CO2",4),rep("CO2+Ca(OH)2",4),rep("CO2+NaOH",4)))
compaired <- list(c("Control", "Ca(OH)2"), 
                  c("Control","CO2"), 
                  c("Control","CO2+Ca(OH)2"),
                  c("Control","CO2+NaOH"))
#fig for Asp 缩写：D,天冬氨酸 酸性，分子量：33.10 
aa_PI_OR=aa_PI[order(aa_PI$AD,decreasing = T),]
D_rich=aa_PI_OR[1:1000,]
D_rich1=aa_PI_OR[1:1000,]
AA_seqs
AA_seqs_TOP_D_NR<-protein_seqs[protein_seqs$Protein %in% rownames(D_rich1), ]
AA_seqs_TOP_D<-AA_seqs[rownames(AA_seqs) %in% rownames(D_rich1), ]
TOP_AD=stri_count_regex(AA_seqs_TOP_D$Sequence, "D")/nchar(AA_seqs_TOP_D$Sequence)*100

data_seq= data.frame(geneid=AA_seqs_TOP_D_NR$Protein,seq= AA_seqs_TOP_D$Sequence,percent=TOP_AD)
data_seq_OR=data_seq[order(data_seq$percent,decreasing = T),]


highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]

# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)


p_highas_colm<-colMeans(p_highas_heatmap)
D_cor <- data.frame(conditionshg,D=p_highas_colm)
D_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Glu缩写：E，谷氨酸 酸性，分子量：147.13 
aa_PI_OR=aa_PI[order(aa_PI$AE,decreasing = T),]
E_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
E_cor <- data.frame(conditionshg,E=p_highas_colm)
E_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG Lys缩写：K，赖氨酸 碱性， 
aa_PI_OR=aa_PI[order(aa_PI$AK,decreasing = T),]
K_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
K_cor <- data.frame(conditionshg,K=p_highas_colm)
K_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG His缩写：H，组氨酸 碱性- 
aa_PI_OR=aa_PI[order(aa_PI$AH,decreasing = T),]
H_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
H_cor <- data.frame(conditionshg,H=p_highas_colm)
H_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG Arg缩写：R，精氨酸 碱性，分子量：174.20 
aa_PI_OR=aa_PI[order(aa_PI$AR,decreasing = T),]
R_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
R_cor <- data.frame(conditionshg,R=p_highas_colm)
R_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#isoleucine 
aa_PI_OR=aa_PI[order(aa_PI$AI,decreasing = T),]
I_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
I_cor <- data.frame(conditionshg,I=p_highas_colm)
I_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#Gly缩写：G，甘氨酸 亲水性，分子量：75.07 
aa_PI_OR=aa_PI[order(aa_PI$AG,decreasing = T),]
G_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
G_cor <- data.frame(conditionshg,G=p_highas_colm)
G_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Ala缩写：A，丙氨酸 疏水性，分子量：89.09 
aa_PI_OR=aa_PI[order(aa_PI$AA,decreasing = T),]
A_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
A_cor <- data.frame(conditionshg,A=p_highas_colm)
A_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#Val缩写：V，缬氨酸 疏水性，分子量：117.15- 
aa_PI_OR=aa_PI[order(aa_PI$AV,decreasing = T),]
V_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
V_cor <- data.frame(conditionshg,N=p_highas_colm)
V_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Leu缩写：L，亮氨酸 疏水性，分子量：131.17- 
aa_PI_OR=aa_PI[order(aa_PI$AL,decreasing = T),]
L_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
L_cor <- data.frame(conditionshg,L=p_highas_colm)
L_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Phe缩写：F，苯丙氨酸 疏水性，分子量：165.19- 
aa_PI_OR=aa_PI[order(aa_PI$AF,decreasing = T),]
F_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
F_cor <- data.frame(conditionshg,F=p_highas_colm)
F_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Trp缩写：W，色氨酸 疏水性，分子量：204.23- 
aa_PI_OR=aa_PI[order(aa_PI$AW,decreasing = T),]
W_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
W_cor <- data.frame(conditionshg,W=p_highas_colm)
W_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Tyr缩写：Y，酪氨酸 亲水性，分子量：181.19- 
aa_PI_OR=aa_PI[order(aa_PI$AY,decreasing = T),]
Y_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
Y_cor <- data.frame(conditionshg,Y=p_highas_colm)
Y_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#TAsn缩写：N，天冬酰胺 亲水性，分子量：132.12-- 
aa_PI_OR=aa_PI[order(aa_PI$AN,decreasing = T),]
N_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
N_cor <- data.frame(conditionshg,N=p_highas_colm)
N_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#Gln缩写：Q，谷氨酰胺 亲水性，分子量：146.15- 
aa_PI_OR=aa_PI[order(aa_PI$AQ,decreasing = T),]
Q_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
Q_cor <- data.frame(conditionshg,Q=p_highas_colm)
Q_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#Met缩写：M，甲硫氨酸 疏水性，分子量：149.21-- 
aa_PI_OR=aa_PI[order(aa_PI$AM,decreasing = T),]
M_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
M_cor <- data.frame(conditionshg,M=p_highas_colm)
M_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Ser缩写：S，丝氨酸 亲水性，分子量：105.095.68
aa_PI_OR=aa_PI[order(aa_PI$AS,decreasing = T),]  
S_rich=aa_PI_OR[1:1000,]          

highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
S_cor <- data.frame(conditionshg,S=p_highas_colm)
S_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#Thr缩写：T，苏氨酸 亲水性，分子量：119.12- 
aa_PI_OR=aa_PI[order(aa_PI$AT,decreasing = T),]     
T_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
T_cor <- data.frame(conditionshg,T=p_highas_colm)
T_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#TCys缩写：C，半胱氨酸 亲水性，分子量：121.16- 
aa_PI_OR=aa_PI[order(aa_PI$AC,decreasing = T),]     
C_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
C_cor <- data.frame(conditionshg,C=p_highas_colm)
C_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


#Pro缩写：P，脯氨酸 疏水性，分子量：115.13 
aa_PI_OR=aa_PI[order(aa_PI$AP,decreasing = T),]     
P_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
P_cor <- data.frame(conditionshg,P=p_highas_colm)
P_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#PI-----
aa_PI_OR=aa_PI[order(aa_PI$PI,decreasing = T),]     
PI=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
PI_cor <- data.frame(conditionshg,PI=p_highas_colm)
PI_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
#lowPI-----
aa_PI_OR=aa_PI[order(aa_PI$PI,decreasing = F),]     
lowPI=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
lowPI_cor <- data.frame(conditionshg,lowPI=p_highas_colm)
lowPI_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
#charge-----
aa_PI_OR=aa_PI[order(aa_PI$charge,decreasing = T),]     
charge=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
charge_cor <- data.frame(conditionshg,charge=p_highas_colm)
charge_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

#negcharge- 
aa_PI_OR=aa_PI[order(aa_PI$charge,decreasing =F),]     
negcharge=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-normalized_hp_d[rownames(normalized_hp_d) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
negcharge_cor <- data.frame(conditionshg,negcharge=p_highas_colm)
negcharge_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

ALL_Amino<- data.frame(condition=conditionshg,D=D_cor[,2],E=E_cor[,2],K=K_cor[,2],H=H_cor[,2],R=R_cor[,2],G=G_cor[,2],A=A_cor[,2],V=V_cor[,2],L=L_cor[,2],I=I_cor[,2],M=M_cor[,2],F=F_cor[,2],W=W_cor[,2],P=P_cor[,2],S=S_cor[,2],T=T_cor[,2],C=C_cor[,2],Y=Y_cor[,2],N=N_cor[,2],Q=Q_cor[,2],PI=PI_cor$PI,lowPI=lowPI_cor$lowPI,charge=charge_cor$charge,negcharge=negcharge_cor$negcharge,CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)



cor_data<- cor(ALL_Amino[,-1])
##计算p值
cor_p <- cor_pmat(cor_data)

##默认绘图square
ggcorrplot(cor_data[1:24,25:29],lab=TRUE,lab_col = "black", lab_size = 2, p.mat = NULL)

cor_amino1<- data.frame(round((cor_data[1:24,25:29]),2))
cor_p1<- data.frame(round(cor_p[1:24,25:29],3))
cor_amino2<- cor_amino1
cor_amino1[abs(cor_amino1)<0.5|cor_p1>0.05]="-"

#fig for  asparitc and serine abundance----
pdf("Fig4 D_abundace.pdf",width=6,height=6) 
ggplot(D_cor,aes(conditionshg,D,color=conditionshg),y="abundance", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
pdf("Fig4 S_abundace.pdf",width=6,height=6) 
ggplot(S_cor,aes(conditionshg,S,color=conditionshg),y="abundance", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#fig for cor between amino and alkalinity----

pdf("Fig4 aminoheatmap100.pdf",width=15,height=3) 
labeledHeatmap(Matrix = t(cor_amino2),
               xLabels = row.names(cor_amino2) ,
               yLabels = colnames(cor_amino2),
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix =t(cor_amino1),
               setStdMargins = FALSE,
               cex.text = 0.5)
dev.off()

# anbundance-------
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

#fig pI distribution -------
pdf("Fig4 mRNA_protein.pdf",width=10,height=6)
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
pdf("Fig4 pi7.5.pdf",width=6,height=6) 
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
pdf("Fig4 pi7-7.5.pdf",width=6,height=6) 
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
pdf("Fig4 pi7.pdf",width=6,height=6) 
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
pdf("Fig4 all.pdf",width=6,height=6) 
ggplot(pi_protein10_m2,aes(conditionshg,mean,color=conditionshg),y="Correlation R2", title =names(traits)[i],palette = "jco", ordered=TRUE)+
  geom_boxplot()+  
  geom_point(size=2)+
  theme_classic()+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level =T,test = t.test)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
#FIG protein enrich------------
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
#bprotein GO-of >7.5----
protein_data2=subset(pi_protein,pi>7.5)

proteinlist2=protein_data2[,1]
ego_green_BP <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "BP", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_GO <-data.frame(ego_green_BP@result) 
proteinlist2=protein_data2[,1]
ego_green_CC <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "CC", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_CC <-data.frame(ego_green_CC@result) 

pdf("p7.5_GO_protein_BP.pdf",width=8,height=3) 
ggplot(data_all_green_GO [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("p7.5_GO_protein_CC.pdf",width=8,height=3) 
ggplot(data_all_green_CC [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#protein kegg of green module-------------
ekp_green2 <- enricher(proteinlist2, 
                       TERM2GENE = pathway2gene, 
                       TERM2NAME = pathway2name, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff =1,
                       pAdjustMethod = "BH",
                       minGSSize =5)
ekp_green_2 <- ekp_green2@result[ekp_green2@result$pvalue <0.05&ekp_green2@result$Count>4,]

#delete pathway about disease
pdf("p7.5_protein_kegg_barplot.pdf",width=10,height=16) 
ggplot(ekp_green_2[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#bprotein GO-of 7-7.5----
protein_data2=subset(pi_protein,pi>7&pi<7.5)
proteinlist2=protein_data2[,1]
ego_green_BP <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "BP", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_GO <-data.frame(ego_green_BP@result) 
proteinlist2=protein_data2[,1]
ego_green_CC <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "CC", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_CC <-data.frame(ego_green_CC@result) 

pdf("p7_7.5_GO_protein_BP.pdf",width=8,height=3) 
ggplot(data_all_green_GO [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("p7_7.5_GO_protein_CC.pdf",width=8,height=3) 
ggplot(data_all_green_CC [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#p7_7.5 kegg of green module-------------
ekp_green2 <- enricher(proteinlist2, 
                       TERM2GENE = pathway2gene, 
                       TERM2NAME = pathway2name, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff =1,
                       pAdjustMethod = "BH",
                       minGSSize =5)
ekp_green_2 <- ekp_green2@result[ekp_green2@result$pvalue <0.05&ekp_green2@result$Count>4,]

#delete pathway about disease
pdf("p7_7.5_protein_kegg_barplot.pdf",width=10,height=16) 
ggplot(ekp_green_2[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#bprotein GO-of <7----
protein_data2=subset(pi_protein,pi<7)
proteinlist2=protein_data2[,1]
ego_green_BP <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "BP", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_GO <-data.frame(ego_green_BP@result) 
proteinlist2=protein_data2[,1]
ego_green_CC <- enrichGO(gene = proteinlist2,
                         keyType = "GID", 
                         OrgDb =Porites.orgdb,
                         ont = "CC", #ALL BP或MF或CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05) 
data_all_green_CC <-data.frame(ego_green_CC@result) 

pdf("p7_GO_protein_BP.pdf",width=8,height=3) 
ggplot(data_all_green_GO [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
pdf("p7_GO_protein_CC.pdf",width=8,height=3) 
ggplot(data_all_green_CC [1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="green"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 
#p7 kegg of green module-------------
ekp_green2 <- enricher(proteinlist2, 
                       TERM2GENE = pathway2gene, 
                       TERM2NAME = pathway2name, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff =1,
                       pAdjustMethod = "BH",
                       minGSSize =5)
ekp_green_2 <- ekp_green2@result[ekp_green2@result$pvalue <0.05&ekp_green2@result$Count>4,]

#delete pathway about disease
pdf("p7_protein_kegg_barplot.pdf",width=10,height=16) 
ggplot(ekp_green_2[1:5,], aes(-log(pvalue), fct_reorder(Description, -log(pvalue)), fill=pvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red',high='red')+
  theme_classic() + ylab(NULL)+ 
  theme(axis.text= element_text(size=9,color="black"),text= element_text(size=9,color="black"))+coord_fixed(ratio=4) 
dev.off() 



save.image("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig4E.RData")


