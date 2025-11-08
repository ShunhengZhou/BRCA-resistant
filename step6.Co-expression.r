### reference this article to construct GRNs. Single-cell transcriptome analysis reveals the dynamics of human immune cells during early fetal skin development

rm(list=ls())
library(Seurat)

load("MCF7_H.rds") ##read the expr matrix
TF_target<-read.table("TF-Target-information.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE) ## hTFTarget 
dim(TF_target)  #1342129
breast_TF_target<-TF_target[grep("breast",TF_target$tissue),] ##only in breast
dim(breast_TF_target) #301734
net<-c();
edge<-c();
exprgene<-rownames(MCF7_H);
for( i in 1:nrow(breast_TF_target)){
	if(length(which(exprgene==breast_TF_target[i,1]))==1 &&
		length(which(exprgene==breast_TF_target[i,2]))==1){  ##with expr########
		if(breast_TF_target[i,1]!=breast_TF_target[i,2]){  ##remove loop########
		    new_edge<-paste(breast_TF_target[i,1],breast_TF_target[i,2],sep=";")
			if(length(intersect(new_edge,edge)) != 1){     ##remove duplicate edge########
				edge<-c(edge,new_edge)                                               ##
				net<-rbind(net,breast_TF_target[i,])
			}	
		}	
	}	
}
dim(net)
write.table(net,"retained_net.txt",sep="\t",quote=FALSE)

####################################################################
Idents(MCF7_H)<-MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids)


corr_method<-"pearson"
net<-read.table("retained_net.txt",sep="\t",stringsAsFactors=FALSE)
calc_correlation<-function(expr_mat){
	# expr_mat<-exp_matrix0d
	result<-c()
	for(i in 1:nrow(net)){  
		TF_exp     <-expr_mat[which(rownames(expr_mat)==net[i,1]),]
		target_exp <-expr_mat[which(rownames(expr_mat)==net[i,2]),]
		if(sum(TF_exp > 0) > 5 && sum(target_exp > 0) > 5){
			tmp_corr<-cor.test(as.numeric(TF_exp),as.numeric(target_exp),method =corr_method)
			result<-rbind(result,c(net[i,1],net[i,2],
			signif(as.numeric(tmp_corr$estimate),3),signif(as.numeric(tmp_corr$p.value),3)) )	
			
		}	
	}
	colnames(result)<-c("TF","target","correlation","p-value")
	return(result)
}


exp_matrixDay0_2<-subset(MCF7_H, idents =c("Day0_2"), invert = FALSE); exp_matrixDay0_2<-as.matrix(exp_matrixDay0_2@assays$RNA@data)
corr_res<-calc_correlation(exp_matrixDay0_2); write.table(corr_res,"Day0-2_net.txt",sep="\t",row.names=FALSE,quote=FALSE)


exp_matrixDay2<-subset(MCF7_H, idents =c("Day2"), invert = FALSE); exp_matrixDay2<-as.matrix(exp_matrixDay2@assays$RNA@data)
corr_res<-calc_correlation(exp_matrixDay2); write.table(corr_res,"Day2_net.txt",sep="\t",row.names=FALSE,quote=FALSE)


exp_matrixDay4<-subset(MCF7_H, idents =c("Day4"), invert = FALSE); exp_matrixDay4<-as.matrix(exp_matrixDay4@assays$RNA@data)
corr_res<-calc_correlation(exp_matrixDay4); write.table(corr_res,"Day4_net.txt",sep="\t",row.names=FALSE,quote=FALSE)


exp_matrixDay7<-subset(MCF7_H, idents =c("Day7"), invert = FALSE); exp_matrixDay7<-as.matrix(exp_matrixDay7@assays$RNA@data)
corr_res<-calc_correlation(exp_matrixDay7); write.table(corr_res,"Day7_net.txt",sep="\t",row.names=FALSE,quote=FALSE)
#############p.adj
#############p.adj
net0d_2<-read.table("Day0-2_net.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net2d<-read.table("Day2_net.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net4d<-read.table("Day4_net.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net7d<-read.table("Day7_net.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
adjustP<-function(net,outfile){
	net$p.adj<-p.adjust(net$p.value,method="bonferroni")          # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",     "fdr", "none"
	sig_corr<-net[intersect(which(net$p.value<0.01), which(net$correlation > 0.2) ),]
	write.table(sig_corr,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
	return(dim(sig_corr))	
}
adjustP(net0d_2,"0d_2_net.p.adj.txt");adjustP(net2d,"2d_net.p.adj.txt");
adjustP(net4d,"4d_net.p.adj.txt");adjustP(net7d,"7d_net.p.adj.txt");

########################degree distribution
rm(list=ls())  ###Day0_2
library(igraph)
network<-read.table("3-0d_2_net.p.adj.txt",header=TRUE,sep="\t",stringsAsFactor=FALSE);
mypi<-cbind(as.character(network[,1]),as.character(network[,2]))
g<-graph.edgelist(mypi)#igraph 
d <- degree(g)
dd = degree.distribution(g, mode = "all", cumulative = FALSE)
degree = 1:max(d)
probability = dd[-1]
nonzero.position = which(probability != 0)
probability = probability[nonzero.position]
degree = degree[nonzero.position]
reg = lm(log(probability) ~ log(degree))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("Fig.S3A.pdf",width=9,height=6)
plot(probability ~ degree, log = "xy", xlab = "Degree of Nodes ", ylab = "Probability", col = "#84a98c" ,main = "Degree Distribution of Day0_2")
curve(power.law.fit,col="#52796f",add=T,lwd=1.4)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(5,0.4,text_infor,cex =1.2)
print(paste("Alpha =", round(alpha, 3)))
print(paste("R square =", round(R.square, 3)))
dev.off()
#######################################################
rm(list=ls())  ###Day2
library(igraph)
network<-read.table("3-2d_net.p.adj.txt",header=TRUE,sep="\t",stringsAsFactor=FALSE);
mypi<-cbind(as.character(network[,1]),as.character(network[,2]))
g<-graph.edgelist(mypi)#igraph 
d <- degree(g)
dd = degree.distribution(g, mode = "all", cumulative = FALSE)
degree = 1:max(d)
probability = dd[-1]
nonzero.position = which(probability != 0)
probability = probability[nonzero.position]
degree = degree[nonzero.position]
reg = lm(log(probability) ~ log(degree))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("Fig.S3B",width=9,height=6)
plot(probability ~ degree, log = "xy", xlab = "Degree of Nodes ", ylab = "Probability", col = "#84a98c" ,main = "Degree Distribution of Day2")
curve(power.law.fit,col="#52796f",add=T,lwd=1.4)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(8,0.3,text_infor,cex =1.2)
print(paste("Alpha =", round(alpha, 3)))
print(paste("R square =", round(R.square, 3)))
dev.off()
#######################################################
rm(list=ls())  ###Day4
library(igraph)
network<-read.table("3-4d_net.p.adj.txt",header=TRUE,sep="\t",stringsAsFactor=FALSE);
mypi<-cbind(as.character(network[,1]),as.character(network[,2]))
g<-graph.edgelist(mypi)#igraph 
d <- degree(g)
dd = degree.distribution(g, mode = "all", cumulative = FALSE)
degree = 1:max(d)
probability = dd[-1]
nonzero.position = which(probability != 0)
probability = probability[nonzero.position]
degree = degree[nonzero.position]
reg = lm(log(probability) ~ log(degree))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("Fig.S3C",width=9,height=6)
plot(probability ~ degree, log = "xy", xlab = "Degree of Nodes ", ylab = "Probability", col = "#84a98c" ,main = "Degree Distribution of Day4")
curve(power.law.fit,col="#52796f",add=T,lwd=1.4)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(20,0.1,text_infor,cex =1.2)
print(paste("Alpha =", round(alpha, 3)))
print(paste("R square =", round(R.square, 3)))
dev.off()
#######################################################
rm(list=ls())  ###Day7
library(igraph)
network<-read.table("3-7d_net.p.adj.txt",header=TRUE,sep="\t",stringsAsFactor=FALSE);
mypi<-cbind(as.character(network[,1]),as.character(network[,2]))
g<-graph.edgelist(mypi)#igraph 
d <- degree(g)
dd = degree.distribution(g, mode = "all", cumulative = FALSE)
degree = 1:max(d)
probability = dd[-1]
nonzero.position = which(probability != 0)
probability = probability[nonzero.position]
degree = degree[nonzero.position]
reg = lm(log(probability) ~ log(degree))
cozf = coef(reg)
power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
pdf("Fig.S3D",width=9,height=6)
plot(probability ~ degree, log = "xy", xlab = "Degree of Nodes ", ylab = "Probability", col = "#84a98c" ,main = "Degree Distribution of Day7")
curve(power.law.fit,col="#52796f",add=T,lwd=1.4)
a = exp(cozf[[1]])
alpha = -cozf[[2]]
R.square = summary(reg)$r.squared
text_infor=paste("R square:",round(R.square,3),sep=" ")
text(10,0.2,text_infor,cex =1.2)
print(paste("Alpha =", round(alpha, 3)))
print(paste("R square =", round(R.square, 3)))
dev.off()

#########Simpson Index


net0d_2<-read.table("0d_2_net.p.adj.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net2d<-read.table("2d_net.p.adj.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net4d<-read.table("4d_net.p.adj.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
net7d<-read.table("7d_net.p.adj.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
shared_nodes<-function(net1,net2){
	node1<-unique(c(net1[,1],net1[,2]))
	node2<-unique(c(net1[,1],net2[,2]))
	simpson<-length(intersect(node1,node2))/min(length(node1),length(node2))
	return(simpson)
	# jaccard<-length(intersect(node1,node2))/length(unique(c(node1,node2)))
	# return(jaccard)
}
shared_edges<-function(net1,net2){
	edge1<-unique(paste(net1[,1],net1[,2],sep=";"))
	edge2<-unique(paste(net2[,1],net2[,2],sep=";"))
	simpson<-length(intersect(edge1,edge2))/min(length(edge1),length(edge2))
	return(simpson)
	# jaccard<-length(intersect(edge1,edge2))/length(unique(c(edge1,edge2)))
	# return(jaccard)
}
result1<-matrix(0,ncol=4,nrow=4)

result1[1,2]<-shared_nodes(net0d_2,net2d);
result1[1,3]<-shared_nodes(net0d_2,net4d);  result1[1,4]<-shared_nodes(net0d_2,net7d);

result1[2,3]<-shared_nodes(net2d,net4d);  result1[2,4]<-shared_nodes(net2d,net7d);

result1[3,4]<-shared_nodes(net4d,net7d);
##########################################
result2<-matrix(0,ncol=4,nrow=4)

result2[1,2]<-shared_edges(net0d_2,net2d);
result2[1,3]<-shared_edges(net0d_2,net4d);  result2[1,4]<-shared_edges(net0d_2,net7d);

result2[2,3]<-shared_edges(net2d,net4d);  result2[2,4]<-shared_edges(net2d,net7d);

result2[3,4]<-shared_edges(net4d,net7d);

result3<-t(result2)

result<-result1+result3 
diag(result)<-1
colnames(result)<-c("Day0_2","Day2","Day4","Day7")
rownames(result)<-c("Day0_2","Day2","Day4","Day7")
library("pheatmap")
pheatmap(result,
         cellheight=30,
         cluster_cols=F,
         cluster_rows=F,
		 filename = "Fig.3E.pdf"	 
)


Dynet_score<-read.csv("DyNet Central Reference Network default node.csv",stringsAsFactors=FALSE,header=TRUE)  ## export from the cytoscape app

Dynet_score<-Dynet_score[order(Dynet_score$DyNet.Rewiring..Dn.score.,decreasing=TRUE),]
Dynet_score<-Dynet_score[1:10,]
data1<-data.frame(Dynet_score$shared.name,Dynet_score$DyNet.Rewiring..Dn.score.)
colnames(data1)<-c("Gene","RewiringScore")

data1 <- within(data1, 
                   Gene <- factor(Gene, levels=data1$Gene  ))

gg<-ggplot(data1,aes(x = Gene, y = RewiringScore))+
	geom_bar(stat="identity",fill="#56A47F")+
	theme_classic()+
	xlab("")+ylab("")+theme( axis.text.y = element_text(size = 16),
   axis.text.x = element_text(angle=90,size = 16,hjust=1,vjust=0.5,colour=rep("#000000",nrow(data1))))
gg


ggsave(gg,file="Fig.3F.pdf")















