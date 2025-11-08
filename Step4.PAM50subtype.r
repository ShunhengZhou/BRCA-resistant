##########PAM50 subtype###
library(genefu)
library(Seurat)
rm(list=ls())

load("MCF7_H.rds") ##read the expr matrix
exp_matrix<-as.matrix(MCF7_H@assays$RNA@data)
annot_data<-read.table("5-sort.txt",header=TRUE,stringsAsFactor=FALSE,sep="\t")  ## ID conversion
##########
common_gene<-intersect(annot_data$probe,rownames(exp_matrix))

exp_matrix<-exp_matrix[common_gene,]
annot_data<-annot_data[annot_data$prob %in% common_gene,]
# dim(exp_matrix)
# dim(annot_data)

gene_exp<-t(exp_matrix)
SubPred_pam50<-molecular.subtyping(sbt.model = "pam50",data = gene_exp,
 annot = annot_data,do.mapping = T)
table(SubPred_pam50$subtype)
write.table(SubPred_pam50$subtype,"7-subtype_PAM50.txt",quote=FALSE,sep="\t")
library(ggplot2)
library(reshape2)  
library(dplyr) 
subtype<-read.table("7-subtype_PAM50.txt",header=TRUE,stringsAsFactor=FALSE,sep="\t")
split_vector <- function(string) strsplit(string, "[.]")
cells<-sapply(rownames(subtype), split_vector)
cells<-matrix(unlist(cells), ncol = 2, byrow = TRUE)

subtype<-cbind(subtype,cells)
colnames(subtype)<-c("subtype","celltype","ID")
freq<-table(subtype$subtype,subtype$celltype)
freq<-prop.table(freq, 2)  ##calculate freq by column
dataset<-setNames(melt(freq), c('Subtype', 'Times', 'percent'))##
write.table(dataset,"7-2subtype_prob.txt",quote=FALSE,sep="\t")
subtype<-read.table("7-2subtype_prob.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
subtype[which(subtype[,2]=="MCF7_RM_CD44H"),2]<-"0d"; subtype[which(subtype[,2]=="MCF7_WM2d_CD44H"),2]<-"2d";
subtype[which(subtype[,2]=="MCF7_WM4d_CD44H"),2]<-"4d"; subtype[which(subtype[,2]=="MCF7_WM7d_CD44H"),2]<-"7d";
##############################
subtype$Subtype <- factor(subtype$Subtype , levels=c("LumA", "LumB","Basal",  "Her2", "Normal") )
subtype$Times <- factor(subtype$Times , levels=c("0d", "2d", "4d", "7d") )
##  ######              
cols<-c("#d0eef0","#a1b2c6","#edd6c7","#b9e5e8","#c7b67a")   
###########################################################		
	
subtype$Time<-rep(c(3,1,5,4,2),each=5)   #######

df_by_type <- group_by(.data = subtype, Time)
df_summarize <- mutate(.data = df_by_type, Fraction = percent/sum(percent))

# df_summarize$Time <- factor(df_summarize$Times , levels=c(3,1,5,4,2) )
gg<-ggplot(data = df_summarize, mapping = aes(x = Time, y = Fraction, fill = Subtype)) + geom_area()  + 
       geom_line(colour = '#2E2E2E', size = 1, position = 'stack', alpha = 0.6) + 
		scale_fill_manual(values=cols)+
		scale_x_continuous(breaks = c(1,2,3, 4,5), labels = c("Day0_1","Day0_2","Day2","Day4","Day7"))+
		geom_vline(xintercept = c(2,3,4),colour="#FFFFFF",linetype="dashed",size=1.6)+
		theme_classic()
		
ggsave(gg,file="Fig.2A.pdf",width=12,height=9)
#






















