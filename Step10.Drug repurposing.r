# ref https://www.jianshu.com/p/218231bba8af  
library(cmapR)
library(Biobase)
library(GeneExpressionSignature) ##PRL algorithm
gct_file<-"L1000/level5_beta_trt_cp_n720216x12328.gctx"
ds_row <- parse_gctx(gct_file, rid=c(1)) ##
MCF7_posi<-grep("MCF7",ds_row@cid)         ##MCF7 data###############
ds_row@cid[MCF7_posi]
colname<-strsplit(ds_row@cid,split = ":")  
small_molecules<-c()
for(i in 1:length(colname)){ 
	if(length(grep("^BRD-",colname[[i]][2]))==1){
		tmp<-strsplit(colname[[i]][2],split = "-")
		small_molecules[i]<-paste(tmp[[1]][1],"-",tmp[[1]][2],sep="") 
	}
	else{
		small_molecules[i]<-colname[[i]][2]
	}
}
write.table(ds_row@cid,file="colnames.txt",sep="\t",quote=FALSE)
unique_SM<-unique(small_molecules)######32869
files=file("L1000/LINCS_small_molecules.tsv",open="r")  ##read by row
LINCS_small_molecules<-c()                             
while ( TRUE ) {
  line = readLines(files, n = 1)
  if ( length(line) == 0 ) {
    break;
  }
  drug<-strsplit(line,split = "\t")
  LINCS_small_molecules<-rbind(LINCS_small_molecules,drug[[1]][1:2])	
}
#################
common_SM<-intersect(LINCS_small_molecules[,1],unique_SM)
# length(unique_SM)  #32869

ds_col <- parse_gctx(gct_file, cid=c(1))  
write.table(ds_col@rid ,file="PRL/rownames.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#########	
for(i in 1:length(common_SM)){   
	posi_SM<-which(small_molecules==common_SM[i])
	posi<-intersect(MCF7_posi,posi_SM)         ##MCF7 cell line##
	posi
	if(length(posi)>0){
		ds <- parse_gctx(gct_file, cid=c(posi))  
			exprs<-ds@mat
		adf<-data.frame(class=c(rep("A",length(posi))));rownames(adf)<-colnames(exprs)
		adf<-new("AnnotatedDataFrame", data = adf)
		exampleSet<-new("ExpressionSet",exprs=exprs,phenoData=adf)
		start=Sys.time()
		MergingSet=RankMerging(exampleSet,"Spearman")   
		end=Sys.time();	end-start

		filename=paste("L1000/PRL/",common_SM[i],sep="")
		write.table(exprs(MergingSet),file=filename,quote=FALSE,row.names=FALSE,col.names=FALSE)
		print(end-start)
		print(length(posi))
	}
}
###############up regulated pattern
library(clusterProfiler)
library(stringr)
library(enrichplot)
library(fgsea)
gsea_sets<-read.csv("mfuzz_3.csv") 
gsea_sets<-gsea_sets[which(gsea_sets[,2]>0.95),]
genes<-gsea_sets[,1]
gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt") ##
allfile<-list.files("L1000/PRL/")

result<-c()
for(i in 1:length(allfile)){
	input_file<-allfile[i]                ##PRL drug response 
	setwd("L1000/PRL/")
	drug_response<-read.table(input_file)
	drug_response<-as.vector(drug_response[,1])

	names(drug_response) <- as.character(gene_annotation$V1)
	geneList <- sort(drug_response, decreasing = TRUE)
	common_gene<-intersect(gene_annotation$V1,gs[,2])
	gsea_sets<-data.frame(ont=rep("MCF7_up_regulated",length(common_gene)),gene=as.numeric(common_gene))
	egmt <- fgsea(pathway=gsea_sets,stats=geneList,scoreType="neg",nperm=10000)
	result<-rbind(result,c(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	print(i)	
}
colnames(result)<-c("drug","ES","NES","pvalue")
write.table(result,file="fgsea_up_regulated.txt",quote=FALSE,sep="\t")

###############down regulated pattern
gsea_sets<-read.csv("mfuzz_6.csv") 
gsea_sets<-gsea_sets[which(gsea_sets[,2]>0.95),]
genes<-gsea_sets[,1]gs <-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_annotation<-read.table("rownames.txt") 
allfile<-list.files()
result<-c()
for(i in 1:length(allfile)){

	input_file<-allfile[i]                ##PRL drug response 
	setwd("L1000/PRL/")
	drug_response<-read.table(input_file)
	drug_response<-as.vector(drug_response[,1])

	names(drug_response) <- as.character(gene_annotation$V1)
	geneList <- sort(drug_response, decreasing = TRUE)
	common_gene<-intersect(gene_annotation$V1,gs[,2])
	gsea_sets<-data.frame(ont=rep("MCF7_down_regulated",length(common_gene)),gene=as.numeric(common_gene))
	class(gsea_sets$gene)
	egmt <- fgsea(pathway=gsea_sets,stats=geneList,scoreType="pos",nperm=10000)
	result<-rbind(result,c(input_file,signif(egmt$ES,3),signif(egmt$NES,3),egmt$pval))
	print(i)	
}
colnames(result)<-c("drug","ES","NES","pvalue")
write.table(result,file="fgsea_down_regulated.txt",quote=FALSE,sep="\t")



###############################WTCS
up<-read.table(file="fgsea_up_regulated-MM-monocle.txt",sep="\t")
colnames(up)<-c("drug_up","ES_up","NES_up","pvalue_up")
down<-read.table(file="fgsea_down_regulated-MM-monocle.txt",sep="\t")
colnames(down)<-c("drug_down","ES_down","NES_down","pvalue_down")

WTCS<-c()
for(i in 1:nrow(up)){
	if(sign(up$NES[i]) == sign(down$NES[i])){
		WTCS[i]<-0
	}
	else{
		WTCS[i]<-(up$NES[i]-down$NES[i])/2
	}
}
merge_result<-cbind(up,down,WTCS)
posi<-intersect(which(merge_result$pvalue_up<0.05),which(merge_result$pvalue_down<0.05))
merge_result<-merge_result[posi,]   ####select the significant###

demo=file("G:/project/4-subpopulation/2-exp_data/4-PRL/LINCS_small_molecules.tsv",open="r")
SM<-c()  
firstline = readLines(demo, n = 1)
while ( TRUE ) {
  line = readLines(demo, n = 1) 
  if ( length(line) == 0 ) {
    break;
  }
  tmp<-unlist(strsplit(line,"\t"))
  SM<-rbind(SM,c(tmp[1:3]))
}
colnames(SM)<-c("ID","pert_name","target")
write.table(SM,file="LINCS_ID_target.txt",sep="\t",row.names=FALSE,quote=FALSE)
drug_info<-c()
for(i in 1:nrow(merge_result)){
	posi<-which(SM[,1]==merge_result$drug_up[i])
	drug_info<-rbind(drug_info,SM[posi,])
}
merge_result<-cbind(merge_result,drug_info)
merge_result <- merge_result[order(merge_result$WTCS),]
write.table(merge_result,"merge_result.WTCS.txt",sep="\t",row.names=FALSE,quote=FALSE)

####bubble plot of WTCS###########################################
library(ggplot2)
data = read.table("E:/Project/4-subpopulation/2-exp_data/4-PRL/L1000/merge_result.WTCS.txt", header=TRUE, sep="\t",stringsAsFactors = F)
data<-data[1:10,]
data<-data[order(data$WTCS,decreasing=FALSE),]
data<-data[1:10,]
data1<-data.frame(data$pert_name,data$WTCS)
colnames(data1)<-c("Drug","WTCS")
data1 <- within(data1,Drug <- factor(Drug, levels=data1$Drug  ))

library(ggplot2)


gg<-   ggplot(data1, aes(x = Drug, y = WTCS)) +ylim(-3,-6)+
  geom_segment(aes(xend = Drug, yend = -3), color = "#84A98C",size = 1.5) +  #
  geom_point(size = 5.5, color = "#84A98C") +  # 
  labs(title = "Lollipop Chart", x = "", y = "") +
  theme_classic()+  
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 18,color="black"),
		axis.title.x = element_text( size = 8),
		axis.text.y = element_text( size = 18,color="black"),
		axis.title.y = element_text( size = 18)
	  )+ 
scale_y_reverse()+ 
 geom_text(aes(label = WTCS),vjust =-0.8,size=5)+
 ylab("WTCS")
ggsave(gg,file="Fig.6B.pdf",width=10,height=6)



#######################barplot of NES
data = read.table("merge_result.WTCS.txt", header=TRUE, sep="\t",stringsAsFactors = F)
data<-data[1:10,]
plt_data = data.frame(
  Class = c(rep("Up_regulated", nrow(data)), rep("Down_regulated", nrow(data))),
  Small_molecules = data$pert_name,
  NES_value = c(data$NES_up, data$NES_down),
  p_value = c(data$pvalue_up, data$pvalue_down),
  stringsAsFactors = F)
plt_data$Small_molecules = factor(plt_data$Small_molecules, levels = c(data$pert_name[1:10])  )
gg<-ggplot(plt_data, aes(x = Small_molecules, y = NES_value, fill = Class)) +
  geom_bar(stat="identity", position="identity",width = 0.48) +scale_fill_manual(values = c( "#84A98C","#F8AD9D"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14,color="black"),
		axis.title.x = element_text( size = 8),
		axis.text.y = element_text( size = 18,color="black"),
		axis.title.y = element_text( size = 18)
	  )+ 
	  ylim(-9,2)+xlab("")+ylab("NES value")+
	geom_text(aes(label = NES_value),vjust =0,size=5)
ggsave(gg,file="Fig.S6A.pdf",width=10,height=6)


######Gene set score
##############################################
library(Seurat)
library(dplyr)
library(ggpubr)
library(reshape2)  
load("MCF7_H.rds")
Idents(MCF7_H)<- MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids) 

########reverse top expressed genes ####################
data = read.table("merge_result.WTCS.txt", header=TRUE, sep="\t",stringsAsFactors = F)
drug<-data$pert_name[2]
drug_ID<-read.table("LINCS_ID_target.txt", header=TRUE, sep="\t",stringsAsFactors = F)
ID<-drug_ID[which(drug_ID$pert_name ==drug),1]
drug_response<-read.table(ID)
drug_response<-as.vector(drug_response[,1])
gene_annotation<-read.table("rownames.txt") ##
names(drug_response) <- as.character(gene_annotation$V1)
drug_response<-sort(drug_response)


##up regulated genes
top_genes<-tail(names(drug_response),200)
gs <-bitr(top_genes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")

MCF7_H <- AddModuleScore(  object = MCF7_H,  features = list(gs$SYMBOL),    name = 'Score')
cols<-c("#516393","#CBCE05","#F07D00","#9D80BA")
plt_data <- data.frame(
	Group <- Idents(MCF7_H),	Score<-as.numeric(MCF7_H$Score1	))
table(plt_data$Group )
colnames(plt_data)<-c("Group","Score")
plt_data <- subset(plt_data , Group = c("Day0_2", "Day2", "Day4", "Day7") )
plt_data$Group <- factor(plt_data$Group , levels=c("Day0_2", "Day2", "Day4", "Day7") )
table(plt_data$Group )
plt_data <-na.omit(plt_data)
table(plt_data$Group )

gg<-ggplot(plt_data, aes(x=Group, y=Score, fill=Group),alpha=0.9) + # fill=name allow to automatically dedicate a color for each group
		geom_violin(draw_quantiles = c( 0.5))+ 
		# geom_violin()+ 
		theme_classic()+
		labs(y = "Score", x="")+
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		stat_compare_means(label.x = 2.4)+
		scale_fill_manual(values=cols)
ggsave(gg,file="Fig.6C.pdf")

## down regulated genes
down_genes<-head(names(drug_response),200)
gs <-bitr(down_genes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
MCF7_H <- AddModuleScore(  object = MCF7_H,  features = list(gs$SYMBOL),    name = 'Score')
cols<-c("#516393","#CBCE05","#F07D00","#9D80BA")
plt_data <- data.frame(
	Group <- Idents(MCF7_H),	Score<-as.numeric(MCF7_H$Score1	))
table(plt_data$Group )
colnames(plt_data)<-c("Group","Score")
plt_data <- subset(plt_data , Group = c("Day0_2", "Day2", "Day4", "Day7") )
plt_data$Group <- factor(plt_data$Group , levels=c("Day0_2", "Day2", "Day4", "Day7") )
table(plt_data$Group )
plt_data <-na.omit(plt_data)
table(plt_data$Group )

gg<-ggplot(plt_data, aes(x=Group, y=Score, fill=Group),alpha=0.9) + # fill=name allow to automatically dedicate a color for each group
		geom_violin(draw_quantiles = c( 0.5))+ 
		theme_classic()+
		labs(y = "Score", x="")+
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		stat_compare_means(label.x = 2.4)+
		scale_fill_manual(values=cols)

ggsave(gg,file="Fig.S6B.pdf")


	












	


















