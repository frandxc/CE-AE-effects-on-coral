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
setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code')
load("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/fig code/RData/Fig2.RData")
list.files("D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data")
AA_seqs <-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/peptide.csv', header = T) 
traits = read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/traits.csv', header = T, row.names=1);
protein_seqs <-read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/protein_LFQ.csv', header = T) 
rownames(AA_seqs)<-AA_seqs $Leading.razor.protein
fa=AA_seqs$Sequence
seqPI=data.frame(fa) %>%
  rownames_to_column("name") %>%
  mutate("length" = Peptides::lengthpep(seq = fa)) %>%  # lengthpep() 计算长度
  mutate("mw" = mw(seq = fa)) %>%                 # mw() 计算分子量
  mutate("hydrophobicity" = hydrophobicity(seq = fa)) %>%     # hydrophobicity() 计算疏水性
  mutate("PI" = pI(seq = fa)) %>%                 # pI() 计算等电点
  mutate("charge" = charge(seq = fa))%>% 
  mutate("instaIndex" = instaIndex(seq = fa))

rownames(seqPI)<-AA_seqs$Leading.razor.protein
DNA_seq=AA_seqs$Sequence
# letterFrequency(DNA_seq[1],"A")
paste(DNA_seq[[1]], collapse="") #extract sequence
nchar(DNA_seq)
#Unigene0168205
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
#fig for Asp 缩写：D,天冬氨酸 酸性，分子量：33.10----
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
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]

# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)


p_highas_colm<-colMeans(p_highas_heatmap)
D_cor <- data.frame(conditionshg,D=p_highas_colm)
D_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Glu缩写：E，谷氨酸 酸性，分子量：147.13----
aa_PI_OR=aa_PI[order(aa_PI$AE,decreasing = T),]
E_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
E_cor <- data.frame(conditionshg,E=p_highas_colm)
E_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG Lys缩写：K，赖氨酸 碱性，------
aa_PI_OR=aa_PI[order(aa_PI$AK,decreasing = T),]
K_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
K_cor <- data.frame(conditionshg,K=p_highas_colm)
K_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG His缩写：H，组氨酸 碱性-----
aa_PI_OR=aa_PI[order(aa_PI$AH,decreasing = T),]
H_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
H_cor <- data.frame(conditionshg,H=p_highas_colm)
H_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG Arg缩写：R，精氨酸 碱性，分子量：174.20-----
aa_PI_OR=aa_PI[order(aa_PI$AR,decreasing = T),]
R_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
R_cor <- data.frame(conditionshg,R=p_highas_colm)
R_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR isoleucine -----
aa_PI_OR=aa_PI[order(aa_PI$AI,decreasing = T),]
I_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
I_cor <- data.frame(conditionshg,I=p_highas_colm)
I_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR Gly缩写：G，甘氨酸 亲水性，分子量：75.07-----
aa_PI_OR=aa_PI[order(aa_PI$AG,decreasing = T),]
G_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
G_cor <- data.frame(conditionshg,G=p_highas_colm)
G_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Ala缩写：A，丙氨酸 疏水性，分子量：89.09-----
aa_PI_OR=aa_PI[order(aa_PI$AA,decreasing = T),]
A_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
A_cor <- data.frame(conditionshg,A=p_highas_colm)
A_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR Val缩写：V，缬氨酸 疏水性，分子量：117.15-----
aa_PI_OR=aa_PI[order(aa_PI$AV,decreasing = T),]
V_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
V_cor <- data.frame(conditionshg,N=p_highas_colm)
V_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Leu缩写：L，亮氨酸 疏水性，分子量：131.17-----
aa_PI_OR=aa_PI[order(aa_PI$AL,decreasing = T),]
L_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
L_cor <- data.frame(conditionshg,L=p_highas_colm)
L_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Phe缩写：F，苯丙氨酸 疏水性，分子量：165.19-----
aa_PI_OR=aa_PI[order(aa_PI$AF,decreasing = T),]
F_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
F_cor <- data.frame(conditionshg,F=p_highas_colm)
F_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Trp缩写：W，色氨酸 疏水性，分子量：204.23-----
aa_PI_OR=aa_PI[order(aa_PI$AW,decreasing = T),]
W_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
W_cor <- data.frame(conditionshg,W=p_highas_colm)
W_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Tyr缩写：Y，酪氨酸 亲水性，分子量：181.19-----
aa_PI_OR=aa_PI[order(aa_PI$AY,decreasing = T),]
Y_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
Y_cor <- data.frame(conditionshg,Y=p_highas_colm)
Y_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR TAsn缩写：N，天冬酰胺 亲水性，分子量：132.12-----
aa_PI_OR=aa_PI[order(aa_PI$AN,decreasing = T),]
N_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
N_cor <- data.frame(conditionshg,N=p_highas_colm)
N_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR Gln缩写：Q，谷氨酰胺 亲水性，分子量：146.15-----
aa_PI_OR=aa_PI[order(aa_PI$AQ,decreasing = T),]
Q_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
Q_cor <- data.frame(conditionshg,Q=p_highas_colm)
Q_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR Met缩写：M，甲硫氨酸 疏水性，分子量：149.21-----
aa_PI_OR=aa_PI[order(aa_PI$AM,decreasing = T),]
M_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
M_cor <- data.frame(conditionshg,M=p_highas_colm)
M_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Ser缩写：S，丝氨酸 亲水性，分子量：105.095.68-----
aa_PI_OR=aa_PI[order(aa_PI$AS,decreasing = T),]  
S_rich=aa_PI_OR[1:1000,]          

highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
S_cor <- data.frame(conditionshg,S=p_highas_colm)
S_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Thr缩写：T，苏氨酸 亲水性，分子量：119.12-----
aa_PI_OR=aa_PI[order(aa_PI$AT,decreasing = T),]     
T_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
T_cor <- data.frame(conditionshg,T=p_highas_colm)
T_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR TCys缩写：C，半胱氨酸 亲水性，分子量：121.16-----
aa_PI_OR=aa_PI[order(aa_PI$AC,decreasing = T),]     
C_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
C_cor <- data.frame(conditionshg,C=p_highas_colm)
C_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)


# FIG FOR Pro缩写：P，脯氨酸 疏水性，分子量：115.13-----
aa_PI_OR=aa_PI[order(aa_PI$AP,decreasing = T),]     
P_rich=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
P_cor <- data.frame(conditionshg,P=p_highas_colm)
P_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR hydrophobicity-----
aa_PI_OR=aa_PI[order(aa_PI$hydrophobicity,decreasing = T),]     
hydrophobicity=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
hydrophobicity_cor <- data.frame(conditionshg,hydrophobicity=p_highas_colm)
hydrophobicity_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
# FIG FOR lowhydrophobicity-----
aa_PI_OR=aa_PI[order(aa_PI$hydrophobicity,decreasing =F),]     
lowhydrophobicity=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
lowhydrophobicity_cor <- data.frame(conditionshg,lowhydrophobicity=p_highas_colm)
lowhydrophobicity_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
# FIG FOR PI-----
aa_PI_OR=aa_PI[order(aa_PI$PI,decreasing = T),]     
PI=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
PI_cor <- data.frame(conditionshg,PI=p_highas_colm)
PI_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
# FIG FOR lowPI-----
aa_PI_OR=aa_PI[order(aa_PI$PI,decreasing = F),]     
lowPI=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
lowPI_cor <- data.frame(conditionshg,lowPI=p_highas_colm)
lowPI_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)
# FIG FOR charge-----
aa_PI_OR=aa_PI[order(aa_PI$charge,decreasing = T),]     
charge=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
charge_cor <- data.frame(conditionshg,charge=p_highas_colm)
charge_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

# FIG FOR negcharge-----
aa_PI_OR=aa_PI[order(aa_PI$charge,decreasing =F),]     
negcharge=aa_PI_OR[1:1000,]
highas=datahgraw1[rownames(aa_PI_OR[1:1000,]),]
p_highas_heatmap<-data_imp_norm[rownames(data_imp_norm) %in% rownames(aa_PI_OR[1:1000,]), ]
highas_P<-highas[rownames(highas) %in% rownames(p_highas_heatmap), ]
# pheatmap(p_highas_heatmap,scale='row',border_color=NA, cluster_cols = FALSE,cluster_rows = T,col=greenred(75),legend =FALSE)
p_highas_colm<-colMeans(p_highas_heatmap)
negcharge_cor <- data.frame(conditionshg,negcharge=p_highas_colm)
negcharge_cor$conditionshg <- factor(conditionshg, level=c("Control", "Ca(OH)2","CO2","CO2+Ca(OH)2","CO2+NaOH"), ordered=TRUE)

ALL_Amino<- data.frame(condition=conditionshg,D=D_cor[,2],E=E_cor[,2],K=K_cor[,2],H=H_cor[,2],R=R_cor[,2],G=G_cor[,2],A=A_cor[,2],V=V_cor[,2],L=L_cor[,2],I=I_cor[,2],M=M_cor[,2],F=F_cor[,2],W=W_cor[,2],P=P_cor[,2],S=S_cor[,2],T=T_cor[,2],C=C_cor[,2],Y=Y_cor[,2],N=N_cor[,2],Q=Q_cor[,2],hydrophobicity=hydrophobicity_cor$hydrophobicity,lowhydrophobicity=lowhydrophobicity_cor$lowhydrophobicity,PI=PI_cor$PI,lowPI=lowPI_cor$lowPI,charge=charge_cor$charge,negcharge=negcharge_cor$negcharge,CO2=traits$CO2,HCO3=traits$HCO3.,CO3=traits$CO32.,pH=traits$pH,alk=traits$Alkalinity)


#fig for cor between amino and alkalinity----

cor_data<- cor(ALL_Amino[,-1])
##计算p值
cor_p <- cor_pmat(cor_data)

##默认绘图square
ggcorrplot(cor_data[1:26,27:31],lab=TRUE,lab_col = "black", lab_size = 2, p.mat = NULL)

cor_amino1<- data.frame(round((cor_data[1:26,27:31]),2))
cor_p1<- data.frame(round(cor_p[1:26,27:31],3))
cor_amino2<- cor_amino1
cor_amino1[abs(cor_amino1)<0.5|cor_p1>0.05]="-"


pdf("aminoheatmap1000.pdf",width=15,height=3) 
labeledHeatmap(Matrix = t(cor_amino2),
               xLabels = row.names(cor_amino2) ,
               yLabels = colnames(cor_amino2),
               colorLabels = FALSE,
               colors = blueWhiteRed (40),
               textMatrix =t(cor_amino1),
               setStdMargins = FALSE,
               cex.text = 0.5)
dev.off()


