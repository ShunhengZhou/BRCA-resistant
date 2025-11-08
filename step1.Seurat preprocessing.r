rm(list=ls())
library(ggplot2)
library(Seurat)
library(dplyr)
library(clustree)
library(org.Hs.eg.db)
library(clusterProfiler)
raw_counts_CD44H<-read.table(file="MCF7_RM_CD44H.No_4981.txt",sep="\t")
raw_counts_CD44H_2d<-read.table(file="MCF7_WM2d_CD44H.No_4999.txt",sep="\t")
raw_counts_CD44H_4d<-read.table(file="MCF7_WM4d_CD44H.No_1649.txt",sep="\t")
raw_counts_CD44H_7d<-read.table(file="MCF7_WM7d_CD44H.No_2129.txt",sep="\t")
############### ID convert###########################
gene.df <- bitr(rownames(raw_counts_CD44H_2d) , fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
gene.df<-gene.df[!duplicated(gene.df[,c('ENSEMBL')]),] 
gene.df<-gene.df[!duplicated(gene.df[,c('SYMBOL')]),] 
get_intersect_gene<-function(expr){
	features <- rownames(expr)[rownames(expr) %in% gene.df$ENSEMBL]   ##common genes##########
	new.gene.df<-gene.df[gene.df$ENSEMBL %in% features,]
	expr <- expr[rownames(expr) %in% features,]
	rownames(expr)<-as.vector(new.gene.df$SYMBOL)
	return(expr)
}
raw_counts_CD44H<-get_intersect_gene(raw_counts_CD44H)
raw_counts_CD44H_2d<-get_intersect_gene(raw_counts_CD44H_2d)
raw_counts_CD44H_4d<-get_intersect_gene(raw_counts_CD44H_4d)
raw_counts_CD44H_7d<-get_intersect_gene(raw_counts_CD44H_7d)

#################
MCF7_H <- CreateSeuratObject(cbind(raw_counts_CD44H, raw_counts_CD44H_2d,raw_counts_CD44H_4d,
         raw_counts_CD44H_7d), project = "MCF7_H", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = MCF7_H@var.genes, npcs = 50, verbose = FALSE) %>% 
	RunUMAP(reduction = "pca", dims = 1:30) %>% 
    RunTSNE(reduction = "pca", dims = 1:30)

MCF7_H@meta.data$stim <- c(rep("0d", ncol(raw_counts_CD44H)), rep("2d", ncol(raw_counts_CD44H_2d)), 
    rep("4d", ncol(raw_counts_CD44H_4d)), rep("7d", ncol(raw_counts_CD44H_7d)))#
save(MCF7_H,file="MCF7_H.rds")
