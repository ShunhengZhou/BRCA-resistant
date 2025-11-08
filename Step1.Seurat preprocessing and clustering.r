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




load("MCF7_H.rds")
head(MCF7_H@meta.data)
pdf("Fig.1B.pdf.pdf",width = 10, height = 10) # color by stim
DimPlot(object = MCF7_H, reduction = "tsne", pt.size = .6, group.by = "stim",cols=c("#409DD8","#CBCE05","#F07D00","#9D80BA"))+
theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
dev.off()
MCF7_H<-FindClusters(MCF7_H,resolution = 0.01)
MCF7_H<-FindClusters(MCF7_H,resolution = 0.05)
#Silhouette score
#reference of this publications https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4
reductions<-"pca"
library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = MCF7_H@reductions$pca))
calc_mean_silhouette<-function(reduction_tag,resolution){
	clusters <- resolution
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	colnames(sil)<-c("cluster", "neighbor", "sil_width")
	# mixing metric
	mat <- data.frame(sil[,])
	mean_silhouette<-mat %>%                      # Specify data frame
	  group_by(cluster) %>%                       # Specify group indicator
	  summarise_at(vars(sil_width),              # Specify column
				   list(name = mean))   
	return(	data.frame(group=rep(reduction_tag,length(mean_silhouette$name)),score=mean_silhouette$name))
}

R0.05<-calc_mean_silhouette("R0.05",MCF7_H$RNA_snn_res.0.05)
R0.1<-calc_mean_silhouette("R0.1",MCF7_H$RNA_snn_res.0.1)
R0.2<-calc_mean_silhouette("R0.2",MCF7_H$RNA_snn_res.0.2)
R0.3<-calc_mean_silhouette("R0.3",MCF7_H$RNA_snn_res.0.3)
R0.4<-calc_mean_silhouette("R0.4",MCF7_H$RNA_snn_res.0.4)
R0.5<-calc_mean_silhouette("R0.5",MCF7_H$RNA_snn_res.0.5)
library(ggplot2)

data<-data.frame(rbind(R0.05,R0.1,R0.2,R0.3,R0.4,R0.5))
colnames(data)<-c("group","score")
cols<-c("#b5e48c","#99d98c","#76c893","#52b69a","#34a0a4","#469d89")

set.seed(1)
gg<-ggplot(data, aes(x=group, y=score)) + 
	labs(y = "Silhouette score", x="Resolution")+
	geom_boxplot(
		color=c("#000000"),
		fill =cols,
		na.rm = TRUE,outlier.shape = NA
	)+
	theme_classic() +
	theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))+ 
	geom_jitter()	
# gg	
ggsave(gg,file='Fig.1C.pdf',width=10,height=10)

pdf("Fig.1D.pdf",width = 10, height = 10) #        
DimPlot(MCF7_H, reduction = "tsne",group.by = 'RNA_snn_res.0.1', label = TRUE, pt.size = .6,
   cols=c("#CBCE05","#409DD8","#9D80BA","#F07D00","#516393"))+
theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
dev.off()


freq<-table(MCF7_H$RNA_snn_res.0.1,MCF7_H$stim) 
freq<-prop.table(freq, 2)  ## 
library(reshape2) 
dataset<-setNames(melt(freq), c('cluster', 'Times', 'percent'))##
write.table(dataset,"1subtype_prob.txt",quote=FALSE,sep="\t")
subtype<-read.table("1subtype_prob.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
######################order of cells########
subtype$cluster <- factor(subtype$cluster , levels=c(0, 1,2, 3, 4) )
subtype$type<-"NA";
subtype$type[which(subtype$cluster==0)]="Day2";
subtype$type[which(subtype$cluster==1)]="Day0_1";
subtype$type[which(subtype$cluster==2)]="Day7";
subtype$type[which(subtype$cluster==3)]="Day4";
subtype$type[which(subtype$cluster==4)]="Day0_2";
subtype$Times <- factor(subtype$Times , levels=c("0d", "2d", "4d", "7d") )
##Statistics of the fraction at each time point######              
cols<-c("#CBCE05","#409DD8","#9D80BA","#F07D00","#516393")   
cols<-c("#409DD8","#516393","#CBCE05","#F07D00","#9D80BA")   
subtype$type <- factor(subtype$type , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )
gg<-ggplot(subtype, aes(fill=type, y=percent, x=Times)) + 
    scale_fill_manual(values=cols)+
	geom_bar( position="fill",stat="identity",width = 0.7)+
	theme_classic()+
	theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(gg,file="Fig.1E.pdf")



MCF7_H.markers <- FindAllMarkers(MCF7_H,  min.pct = 0.25, logfc.threshold = 0.25)
top10markers<-MCF7_H.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)  ##top 10 of each cluster 
top10<-top10markers$gene[c(11:20,41:50,1:10,31:40,21:30)]
markergene_list<-unique(as.character(top10))
pdf("Fig.S1.pdf",width = 10, height = 10) 
DoHeatmap(MCF7_H, 
    angle = 45,size=8, features = markergene_list,group.colors=c("#409DD8","#516393","#CBCE05","#F07D00","#9D80BA")  
   )+scale_fill_gradientn(colors = c("#50114F", "#000000", "#FDFD0F"))+ theme(axis.text.y = element_text(size = 10))
dev.off()  


	