rm(list=ls())
install.packages("oncoPredict")
install.packages("Seurat")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
#https://cran.r-project.org/web/packages/oncoPredict/vignettes/calcPhenotype.html
library(oncoPredict)
library(Seurat)
options(stringsAsFactors = F)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(ggplot2)
# 
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

GDSC2_Expr = readRDS('GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res  = readRDS("GDSC2_Res.rds")
GDSC2_Res  <- exp(GDSC2_Res) 
dim(GDSC2_Res)

library(limma)
drugs<-colnames(GDSC2_Res);
drugs<-strsplit2(drugs,split= "_")[,1:2]
drugs[,1]<-toupper(as.character(drugs[,1]))
AI<-c("ANASTRAZOLUM","ANASTROZOLE","ARIMIDEX","ANASTROZOLUM","ANASTRAZOLE","","EXEMESTANE","AROMATASE ","AROMASIN","","","LETROZOLE","FEMARA","","","TAMOXIFEN","TAMOXIPHENE","TAMOXIFEN","NOLVADEX","NOVADEX")
drug<-intersect(AI,drugs[,1]);posi<-which(drugs[,1]==drug)
GDSC2_Res<-GDSC2_Res[,c(1,posi)] ##only Tamoxifen is included in this dataset

load("MCF7_H.rds") ##read the expr matrix
Idents(MCF7_H)<-MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids)
exp_matrix<-as.matrix(MCF7_H@assays$RNA@data)
testExpr<-exp_matrix

dim(GDSC2_Res)
dim(testExpr) 
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


cols<-c("#409DD8","#516393","#CBCE05","#F07D00","#9D80BA") #

IC50<-read.csv("calcPhenotype_Output/DrugPredictions.csv",row.names=1,header=TRUE)
# IC50$Tamoxifen_1199	

data <- data.frame(
	group <- Idents(MCF7_H),	IC50<-IC50$Tamoxifen_1199	)
colnames(data)<-c("group","IC50")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )
gg <- ggplot(data, aes(x=group, y=log(IC50), fill=group),alpha=0.9) + # fill=name allow to automatically dedicate a color for each group
		geom_violin(draw_quantiles = c( 0.5))+ 
		theme_classic()+
		labs(y = "logIC50", x="Time")+
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		stat_compare_means(label.x = 1.4)+
		scale_fill_manual(values=cols)

ggsave(gg,file="Fig.2G.pdf")
			  
			  
			  
	