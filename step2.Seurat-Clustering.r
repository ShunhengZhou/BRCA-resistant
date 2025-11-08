rm(list=ls())
library(Seurat)
library(dplyr)
load("MCF7_H.rds")
head(MCF7_H@meta.data)
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
ggsave(gg,file='4-Silhouette score.pdf',width=10,height=10)










	