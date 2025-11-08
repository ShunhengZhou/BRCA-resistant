# https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# https://cole-trapnell-lab.github.io/monocle3/docs/starting/#cell_data_set
# http://cole-trapnell-lab.github.io/monocle-release/monocle3/   
# https://www.jianshu.com/p/c402b6588e17
rm(list=ls())
library(monocle3)
library(Seurat)
library(ggplot2)
require(dplyr)
load("MCF7_H.rds")
Idents(MCF7_H)<-MCF7_H$RNA_snn_res.0.1
new.cluster.ids <- c("Day2","Day0_1","Day7","Day4","Day0_2") ##set the cell type of each cluster,merge cluster
names(new.cluster.ids) <- levels(MCF7_H)
MCF7_H <- RenameIdents(MCF7_H, new.cluster.ids)
##remove Day0_1
MCF7_H<-subset(MCF7_H,idents=c("Day0_2","Day2","Day4","Day7"))

expression_matrix <- as(as.matrix(MCF7_H@assays$RNA@counts), 'sparseMatrix')
cell.type<-Idents(MCF7_H)
cell_metadata <- cbind(MCF7_H@meta.data,cell.type)
gene_annotation<-data.frame(gene_short_name = row.names(MCF7_H), row.names = row.names(MCF7_H))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method="UMAP")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(MCF7_H, reduction = "tsne")  ##
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
cds <- learn_graph(cds,close_loop = FALSE)

######################
plot_cells(cds, color_cells_by = "cell.type", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)

##############################################
pdf("Fig.5A.pdf")
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()

######################### monocle3 DEG###############
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
write.table(genes,"genes.txt",sep="\t",quote=FALSE)

genes<-read.table("genes.txt",sep="\t")
genes<-as.vector(genes[,1])
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
labels<-c("E2F1","FOXM1","FOS","HSPA5","CD44","NQO1","BCAS3","VMP1","TPX2","CKS1B","PRC1")  

critical_genes<-c("E2F1","FOXM1","BRCA1","TFAP2C","RAD21")  ##top10 和 down-regulated模式交集的5个基因
labels<-union(labels,critical_genes)


labels<-intersect(labels,genes)
posi<-match(labels,genes)


lab = rowAnnotation(ano = anno_mark(at = posi,
labels = labels,
labels_gp = gpar(fontsize = 12)))

htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,show_row_dend = FALSE,
  show_column_names            = FALSE,
  # row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,right_annotation = lab)

pdf("Fig.2B.pdf")
print(htkm)
dev.off()


  
  
  