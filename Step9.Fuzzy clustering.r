##R  3.6.0
#REF https://mp.weixin.qq.com/s/ubfTKFyUeBR585EfR9nfig

library("Mfuzz")
library("Seurat")

load("MCF7_H.rds") ##read the expr matrix
Idents(MCF7_H)<-MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids)
##remove Day0_1##############
MCF7_H<-subset(MCF7_H,idents=c("Day0_2","Day2","Day4","Day7"))
MCF7_H@active.ident <- factor(MCF7_H@active.ident, 
                            levels=c("Day0_2","Day2","Day4","Day7") )

exp_matrix<-as.matrix(MCF7_H@assays$RNA@data)

label<-Idents(MCF7_H)  
df1<-data.frame(label,t(exp_matrix))
##average###
df2<-aggregate(df1[,colnames(df1)[2:ncol(df1)]],by=list(df1$label),mean,na.rm= TRUE)

row.names(df2)<-df2[,1]
DEGs_exp_averp<-data.frame(t(df2[,-1]))
df3a<-as.matrix(DEGs_exp_averp)
df3Ex<- ExpressionSet(assayData = df3a)
df3F <- filter.NA(df3Ex,thres = 0.25)
df3F <- fill.NA(df3F,mode = 'mean')
df3F <- filter.std(df3F,min.std=0.1)

df3F <- standardise(df3F)
set.seed(2022)
cluster=6
cl <- mfuzz(df3F,c=cluster,m=1.25)

cl$cluster[cl$cluster == 1] # 


## visualization
library(RColorBrewer)  
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
mycol <- c("cyan","yellow","orangered"); color.2 <- colorRampPalette(mycol)(100)
write.csv(cl$membership,file="membership.csv") 
dir.create(path="2mfuzz",recursive = TRUE)
for(i in 1:c){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("2mfuzz","/mfuzz_",i,".csv"))
}
############DIY###################
color<-c("#a4c3b2","#C9AFD4","#f8ad9d","#c9cba3","#83c5be","#84a98c")  ##
for(i in 1:6){
	cluster_No<-i;
	cluster_color<-color[i]
	library(reshape2)
	cluster_gene<-intersect(names(which(cl$cluster==cluster_No)),names(which(cl$membership[,cluster_No]>0.95)))
	cluster_expr<-df3F@assayData$exprs[cluster_gene,]
	#########average
	mean_expr<-colMeans(cluster_expr);
	mean_dataset<-data.frame(gene="mean",Days=c("Day0_2","Day2","Day4","Day7"),Expression=mean_expr,label="mean")

	dataset<-setNames(melt(cluster_expr), c('gene', 'Days', 'Expression'))
	dataset$label<-"cluster"
	dataset<-rbind(dataset,mean_dataset)
	gg<-ggplot(dataset, aes(x = Days, y = Expression, color = label,group=gene)) +
		geom_line(aes(size = label))+
		scale_color_manual(values=c(cluster_color, "black"))+
		scale_size_manual(values = c(cluster = 0.618, mean = 1.5))+
		theme_classic()+
		theme(axis.text.y = element_text(size = 25),axis.text.x = element_text(size = 28),axis.title = element_text(size = 30))+
		theme(legend.position="none")
	filename<-paste("Fig.5C.","cluster",cluster_No,".pdf",sep="")
	ggsave(gg,file=filename,height=4,width=12)
}


#################enrichment analysis 
library(clusterProfiler)
library(dplyr)
library('org.Hs.eg.db')
library(ggplot2)
clusters<-read.csv(file="mfuzz_3.csv",header=TRUE,stringsAsFactors=FALSE)
setwd("mfuzzenrichment/")
DEgenelist<- clusters[,1]
gene.df <- bitr(DEgenelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
eKEGG<- enrichKEGG(gene.df$ENTREZID, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.05)
pdf(file="Fig.S5A.pdf",width=8,height=8)
p2<-dotplot(eKEGG, showCategory = 10,title="Top10 enriched KEGG pathway")
print(p2)
dev.off() 

clusters<-read.csv(file="mfuzz_6.csv",header=TRUE,stringsAsFactors=FALSE)
DEgenelist<- clusters[,1]
gene.df <- bitr(DEgenelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
eKEGG<- enrichKEGG(gene.df$ENTREZID, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.05)
pdf(file="Fig.S5B.pdf",width=8,height=8)
p2<-dotplot(eKEGG, showCategory = 10,title="Top10 enriched KEGG pathway")
print(p2)
dev.off() 

######################## heatmap of genes in cluster3, 6
load("MCF7_H.rds")
Idents(MCF7_H)<-MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids)  #combine and rename Day0_1 clusters 
MCF7_H<-subset(MCF7_H,idents=c("Day0_2","Day2","Day4","Day7"))

clust3<-read.csv(file="mfuzz_3.csv",header=TRUE,row.names=1)
clust6<-read.csv(file="mfuzz_6.csv",header=TRUE,row.names=1)
markergene_list<-as.character(c(rownames(clust3),rownames(clust6)))
all.genes<- rownames(GetAssayData(object = MCF7_H, slot = "scale.data"))
all.genes<-intersect(markergene_list,all.genes)
write.table(all.genes,"all.genes.txt",row.names=FALSE,quote=FALSE,sep="\t")
genes.to.label<-c("SQSTM1","CDKN1A", "SMAD3","SFN","CDC20","CCNA2","PTTG1","CDK1","CCNB2","PLK1")
critical_genes<-c("E2F1","FOXM1","BRCA1","TFAP2C","RAD21")  ##
genes.to.label<-union(genes.to.label,critical_genes)
labels <- rep(x = "transparent", times = length(x = all.genes))
labels[match(x = genes.to.label, table = all.genes)] <- "black"

all.genes[which(all.genes=="TFAP2C" )]<-all.genes[length(all.genes)]; all.genes[length(all.genes)]<-"TFAP2C" ;            #move TFAP2C to the end, or covered by other gene###################

MCF7_H@active.ident <- factor(MCF7_H@active.ident, levels=c("Day0_2","Day2","Day4","Day7"))

pdf("Fig.5D.pdf",width = 10, height = 10) 

DoHeatmap(MCF7_H, 
    angle = 45,size=6, features = all.genes,slot = "scale.data",,group.colors=c("#516393","#CBCE05","#F07D00","#9D80BA")
   )+scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(axis.text.y = element_text(size = 16,color = rev(x = labels)))
dev.off()  
