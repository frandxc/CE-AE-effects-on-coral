# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c('ggdendro'))

# load package------------
# load library------------
library(AnnotationHub)
library(biomaRt)
library(clusterProfiler)
library("topGO")
library("Rgraphviz")
library("pathview")
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

setwd('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/make orgdb_sym')
anno_all <- read.csv('D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/All_Unigene.advance.annotation.csv', header = T, row.names=1)
# sym blast results----
file_path <- "D:/my paper/2018-OLIVINE/滨共表达和论文/2021R - 副本/raw data/symbiont_result/sym"
files <- list.files(path = file_path , full.names = TRUE)
data_sym <- read.csv(paste0(file_path, "/clade_A1.blast.csv"), sep = "\t", header = T)
# 循环读取剩余两个CSV文件并进行rbind拼接
for (i in 2:12) {
  file <- read.csv(files[i], sep = "\t", header = T)
  data_sym <- rbind(data_sym, file)
}
data_sym_filter <- data_sym[data_sym$pident>90,]
# 假设你的数据框名为data，列名为col_name
data_sym_filter$qaccver <- gsub("\\.1$", "", data_sym_filter$qaccver)
row.names(data_sym_filter) <- make.unique(data_sym_filter$qaccver)
anno_all_sym <- anno_all[row.names(anno_all) %in% row.names(data_sym_filter),]

egg <- anno_all_sym
egg[egg=="--"]<-NA #这个代码来自花花的指导(将空行变成NA，方便下面的去除)
egg[egg=="-"]<-NA

#make orgdb: gene_info----------------
egg_1<-egg[,1:3]
egg_1[,1] <-rownames(egg)
egg_1[,2:3] <- egg[,c(1,5)]
rownames(egg_1) <- NULL
colnames(egg_1) <-c('GID','symbol','nrid')
gene_info <- egg_1 %>%
  dplyr::select(GID =GID,GENENAME =nrid) %>% na.omit()



#make orgdb: gene2go_1----------------
gobind<-egg[,1:2]
gobind[,1] <-rownames(egg)
gobind[,2:4]<- egg[,c(23:25)]
rownames(gobind) <- NULL
colnames(gobind) <-c('GID','GO1','GO2','GO3')
gobind_1 <-data.frame(paste(gobind[,2],gobind[,3],gobind[,4]))

gobind_2 <-gsub("NA","",gobind_1[,1])
gobind_2 <-gsub(" ","",data.frame(gobind_2)[,1])
gobind_2 <-data.frame(gobind_2)

gobind_3<-egg[,1:2]
gobind_3[,1] <-rownames(egg)
gobind_3[,2]<- gobind_2[,1]
colnames(gobind_3) <-c('GID','GO')
gobind_3[gobind_3==""]<-NA

gterms <- gobind_3%>%dplyr::select(GID =GID,GO=GO) %>% na.omit()
gene2go <- data.frame(GID = character(),
                      GO = character())

all_go_list = str_split(gterms$GO,";")
gene2go <- data.frame(GID = rep(gterms$GID,
                                times = sapply(all_go_list,length)),
                      GO = unlist(all_go_list),
                      EVIDENCE = "IEA")

all_go_list_1 = str_split_fixed(gene2go$GO,"//",2)
all_go_list_1 = data.frame(all_go_list_1)
gene2go_1 <- data.frame(GID =rep(gterms$GID,
                                 times = sapply(all_go_list,length)),
                        GO = all_go_list_1$X1,
                        EVIDENCE = "IEA")



#make orgdb: kegg--------------
colnames(egg)
egg_3<-egg[,1:2]
egg_3[,1] <-rownames(egg)
egg_3[,2] <- egg[,c(18)]
rownames(egg_3) <- NULL
names(egg)
colnames(egg_3) <-c('GID','KO')
gene2ko <- egg_3 %>% dplyr::select(GID =GID, KO =KO) %>%na.omit()


all_KO_list = str_split_fixed(gene2ko$KO,"//",2)
all_KO_list = data.frame(all_KO_list)
gene2KO_1 <- data.frame(GID = gene2ko$GID,
                        KO = all_KO_list$X1)

#make orgdb: pathway2name, ko2pathway-----
options(stringsAsFactors = F)
update_kegg <- function(json = "smin00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
}

update_kegg(json = "smin00001.json")

load(file = "kegg_info.RData")

colnames(ko2pathway)=c("KO",'Pathway')

gene2pathway <- gene2KO_1 %>% left_join(ko2pathway,by = "KO") %>% 
  dplyr::select(GID,Pathway) %>% na.omit()


# 查询物种的Taxonomy，例如要查sesame
# https://www.ncbi.nlm.nih.gov/taxonomy/?term=sesame
tax_id = "627025"
genus = "Porites" 
species = "_pukoensis"
library(dplyr)
gene_info <- dplyr::distinct(gene_info)
gene2GO <- dplyr::distinct(gene2go_1)
gene2KO <- dplyr::distinct(gene2KO_1)


ggO <- makeOrgPackage(gene_info=gene_info,
               go=gene2GO,
               ko=gene2KO,
               maintainer='yuanxc<77237211@qq.com>',
               author='yuanxc<77237211@qq.com>',
               pathway=gene2pathway,
               version="0.1",
               outputDir = ".",
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")
head(gene2GO)
