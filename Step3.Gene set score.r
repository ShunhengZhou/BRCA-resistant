rm(list=ls())
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggpubr)
library(reshape2)  
load("MCF7_H.rds")
Idents(MCF7_H)<- MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids) 
MCF7_H$type<-Idents(MCF7_H)


#########cell cycle##############################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
MCF7_H<- CellCycleScoring(MCF7_H, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pdf("Fig.2C.pdf",width = 10, height = 10)
DimPlot(MCF7_H, reduction = "tsne",group.by = 'Phase', label = TRUE, pt.size = .6,cols=c("#C5BFCA","#F0BB41","#8DE15B"))
dev.off()
################
MCF7_H <- AddModuleScore(  object = MCF7_H,  features = list(c(s.genes,g2m.genes)),    name = 'cyclingScore')
p1 <- FeaturePlot(MCF7_H, reduction="tsne",pt.size = 1,features = c("cyclingScore1"), combine = FALSE )
fix.sc <- scale_color_gradientn( colours = c("#185695","#FFFFFF", "#861627"),limits = c(-0.75,0.75))   

p2 <- lapply(p1, function (x) x + fix.sc)
p2 <- lapply(p2, function(x){x + labs(title = "cycling score")})
pdf("Fig.2D.pdf",width=9,height=9)
CombinePlots(p2,ncol=1)
dev.off()


########EMT gene set############## download from dbemt database
cols<-c("#409DD8","#516393","#CBCE05","#F07D00","#9D80BA") #c("#CBCE05","#409DD8","#9D80BA","#F07D00","#516393","#516393")
EMT<-read.table(file="dbemt2.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
MCF7_H <- AddModuleScore(  object = MCF7_H,  features = list(EMT$GeneSymbol),    name = 'EMT')

data <- data.frame(
	group <- Idents(MCF7_H),	score<-MCF7_H$EMT1)
colnames(data)<-c("group","EMT score")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )
gg <- ggplot(data, aes(x=group, y=score, fill=group),alpha=0.9) + # fill=name allow to automatically dedicate a color for each group
		geom_violin()+ 
		theme_classic()+
		labs(y = "EMT score", x="Time")+
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		stat_compare_means(label.x = 1.4)+
		scale_fill_manual(values=cols)
ggsave(gg,file="Fig.S2A.pdf")

##########cell cycle phase ###########

subtype<-data.frame(MCF7_H$Phase,MCF7_H$type)
freq<-table(subtype$MCF7_H.Phase,subtype$MCF7_H.type)
freq<-prop.table(freq, 2)  ##frequence of each subtype###   
dataset<-setNames(melt(freq), c('Phase', 'Times', 'percent'))##
write.table(dataset,"6-2subtype_prob.txt",quote=FALSE,sep="\t")
subtype<-read.table("6-2subtype_prob.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")
######################order of cell cycle########
subtype$Subtype <- factor(subtype$Phase , levels=c("G1", "S","G2M") )
subtype$Times <- factor(subtype$Times , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )
###########################################################		
subtype$Time<-rep(c(3,1,5,4,2),each=3)   ##

df_by_type <- group_by(.data = subtype, Time)
df_summarize <- mutate(.data = df_by_type, Fraction = percent/sum(percent))
gg<-ggplot(data = df_summarize, mapping = aes(x = Time, y = Fraction, fill = Subtype)) + geom_area()  + 
       geom_line(colour = '#2E2E2E', size = 1, position = 'stack', alpha = 0.6) + 
		scale_fill_manual(values=c("#C5BFCA","#F0BB41","#8DE15B"))+
		scale_x_continuous(breaks = c(1,2,3,4,5), labels = c("Day0_1","Day0_2", "Day2", "Day4", "Day7"))+
		geom_vline(xintercept = c(2,3,4),colour="#FFFFFF",linetype="dashed",size=1.6)+
		theme_classic()
		
ggsave(gg,file="Fig.2B.pdf",width=12,height=9)

################################################
data <- data.frame(
	group <- MCF7_H$type,	score<-MCF7_H$cyclingScore1)
colnames(data)<-c("group","cellcycle")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )


gg <- ggplot(data, aes(x=group, y=cellcycle, fill=group)) + # fill=name allow to automatically dedicate a color for each group
		geom_violin()+ 
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		labs(y = "cell cycle score", x="Time")+
		stat_compare_means(label.x = 1.4)+
		scale_fill_manual(values=cols)+
		theme_classic()		
ggsave(gg,file="Fig.2E.pdf")

#############hsa04210 apoptosis pathway  ###
library(KEGGREST)
query <- keggGet(c("hsa04210"))  
names(query[[1]])
genelist<-unlist(lapply(query[[1]]$GENE,function(x) strsplit(x,";")))
apoptosis<-genelist[1:length(genelist)%%3 ==2]

MCF7_H <- AddModuleScore(object = MCF7_H,  features = list(apoptosis),    name = 'apoptosis')
data <- data.frame(
	group <- MCF7_H$type,	score<-MCF7_H$apoptosis1)
colnames(data)<-c("group","score")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )


gg<-ggplot(data, aes(x=group, y=score)) + 
	labs(y = "apoptosis score", x="Time")+
	geom_boxplot(
		color=c("#C0C0C0"),
		fill =cols,
		na.rm = TRUE
	)+
	stat_compare_means(label.x = 1.4)+
	theme_classic() 
ggsave(gg,file="Fig.S2B.pdf")

########stemness pathway### Signaling pathways regulating pluripotency of stem cells 
library(KEGGREST)
query <- keggGet(c("hsa04550"))  ## 
names(query[[1]])
genelist<-unlist(lapply(query[[1]]$GENE,function(x) strsplit(x,";")))
stemness<-genelist[1:length(genelist)%%3 ==2]

MCF7_H <- AddModuleScore(object = MCF7_H,  features = list(stemness),    name = 'stemness')
data <- data.frame(
	group <- MCF7_H$type,	score<-MCF7_H$stemness1)
colnames(data)<-c("group","stemness")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )

gg <- ggplot(data, aes(x=group, y=stemness, fill=group)) + # fill=name allow to automatically dedicate a color for each group
		geom_violin()+ 
		theme_classic()+
		theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +      # 45 
		labs(y = "stemness score", x="Time")+
		stat_compare_means(label.x = 1.4)+
		scale_fill_manual(values=cols)
ggsave(gg,file="Fig.S2C.pdf")

#################ferroptosis gene############################
cols<-c("#409DD8","#516393","#CBCE05","#F07D00","#9D80BA")
ferroptosis<-read.table(file="Ferrdball.txt",sep="\t",stringsAsFactors=FALSE)
suppressor<-ferroptosis[grep("suppressor",ferroptosis$V2),] ##only in breast
MCF7_H <- AddModuleScore(  object = MCF7_H,  features = list(suppressor$V1),    name = 'ferroptosis_suppressor')

# suppressor ################################### 
data <- data.frame(
	group <- MCF7_H$type,	score<-MCF7_H$ferroptosis_suppressor1)
colnames(data)<-c("group","ferroptosis_suppressor scor")
data$group <- factor(data$group , levels=c("Day0_1","Day0_2", "Day2", "Day4", "Day7") )

gg<- ggplot(data, aes(x=group, y=score, fill=group)) + 
	geom_boxplot(fill =cols	)+
	labs(y = "ferroptosis suppressor score", x="Time")+
	stat_compare_means(label.x = 1.4)+
	theme_classic()		
		
ggsave(gg,file="Fig.S2D.pdf")
####################



