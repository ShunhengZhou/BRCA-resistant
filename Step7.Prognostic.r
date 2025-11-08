######################################################################################
#4.1.3  https://stackoverflow.com/questions/34349467/php-variable-in-r
rm(list=ls())
library(survival)
library(survminer)
genes<-c("E2F1","CEBPB","BRCA1","TFAP2C","JUN","GATA3","EGR1","NR2F2","RAC3","C21orf58")
###########TCGA RNA-seq of BRCA###########
expr     <-read.table("TCGA-BRCA.htseq_fpkm.tsv",stringsAsFactors=FALSE,sep="\t")
colnames(expr)<-expr[1,];expr<-expr[-1,]
rownames(expr)<-expr[,1];expr<-expr[,-1]
Tumor_sample<-grep(".01A$",colnames(expr),value=TRUE)
Tumor_expr<-expr[,colnames(expr) %in% Tumor_sample]  ###select the tumor samples##
##########annotation####	
anno_data<-read.table("gencode.v22.annotation.gene.probeMap",header=TRUE,sep="\t",stringsAsFactors=FALSE)
# genes<-c("E2F1","CEBPB","BRCA1","TFAP2C","JUN","GATA3","EGR1")
anno_data<-anno_data[anno_data$gene %in% genes,] 
rownames(anno_data)<-anno_data$id

# Comprehensive molecular portraits of human breast tumors  table S1 to get the subtype information, get the sample barcode of LumA
BRCAsubtype<-read.table("BRCA_subtype.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
BRCAsubtype<-BRCAsubtype[BRCAsubtype$PAM50.mRNA %in% c("Luminal A"),]
barcode<-BRCAsubtype$Complete.TCGA.ID
#########
colname<-gsub("-01A$", "",colnames(Tumor_expr)) 
colnames(Tumor_expr)<-colname
# LumA_expr<-Tumor_expr
LumA_expr<-Tumor_expr[,colname %in% barcode]  ##LumA samples
clinical<-read.table("survival_BRCA_survival.txt",stringsAsFactors=FALSE,sep="\t",fill=TRUE)
colnames(clinical)<-clinical[1,];	clinical<-clinical[-1,]
clinical<-clinical[clinical[,1] %in% grep(".01$",clinical[,1],value=TRUE), ]
rownames(clinical)=clinical[,2]
##############################id convert
LumA_expr<-LumA_expr[rownames(LumA_expr) %in% anno_data$id,]
geneSymbol<-c()
for(i in 1:nrow(LumA_expr)){
	geneSymbol[i]<-anno_data[which(anno_data$id == rownames(LumA_expr)[i]),2]
}
rownames(LumA_expr)<-geneSymbol;
################
common_samples<-intersect(rownames(clinical),colnames(LumA_expr))
a=0;b=0;
for(i in 1:length(common_samples)){   ##sort ###
	a[i]=which(rownames(clinical)==common_samples[i])  ##clinData
	b[i]=which(colnames(LumA_expr)==common_samples[i])##profileData
}
clinical<-clinical[a,]
LumA_expr<-LumA_expr[,b]              ##
####################
getPvalue<-function(x=Figer,digits = max(options()$digits - 4, 3)){### get p value

   if (is.matrix(x$obs)){
	    otmp <- apply(x$obs,1,sum)
	    etmp <- apply(x$exp,1,sum)
    }else{
	    otmp <- x$obs
	    etmp <- x$exp
	    }
   df <- (sum(1*(etmp>0))) -1
   pvalue<-format(signif(1-pchisq(x$chisq, df),digits))
  return(pvalue);
}
#####
#################################cox regression
test4 <- list(time=as.numeric(clinical$OS.time),  
              status=as.numeric(clinical$OS), 
              x1=as.numeric(LumA_expr[1,]),
              x2=as.numeric(LumA_expr[2,]),
              x3=as.numeric(LumA_expr[3,]),
			  x4=as.numeric(LumA_expr[4,]),
              x5=as.numeric(LumA_expr[5,]),
			  x6=as.numeric(LumA_expr[6,]),
              x7=as.numeric(LumA_expr[7,]),
              x7=as.numeric(LumA_expr[8,]),
              x7=as.numeric(LumA_expr[9,]),
              x7=as.numeric(LumA_expr[10,])
) 
t4<-coxph(Surv(time, status) ~ x1+x2+x3+x4+x5+x6+x7, test4) 
rr<-c(summary(t4)[[7]][1],
      summary(t4)[[7]][2],
      summary(t4)[[7]][3],
      summary(t4)[[7]][4],
      summary(t4)[[7]][5],
      summary(t4)[[7]][6],
      summary(t4)[[7]][7], 
      summary(t4)[[7]][8], 
      summary(t4)[[7]][9], 
      summary(t4)[[7]][10] 
	  )
xx<-c()
for (i in 1:ncol(LumA_expr)){
	xx_temp<-0;
	for(j in 1:length(rr)){
		xx_temp<-xx_temp+rr[j]*as.numeric(LumA_expr[j,i])
	}
	xx<-c(xx,xx_temp)
}
gene_order<-rownames(LumA_expr)
write.table(cbind(gene_order,rr),"cox.coefficient.txt",sep="\t",quote=FALSE)
median<-median(as.numeric(xx))
clas<-ifelse(as.numeric(xx) > median ,"High-risk","Low-risk")
survdata<-cbind(clinical,class=clas)


###############cutoff  #########
merge_data<-data.frame(time=as.numeric(survdata$OS.time),event=as.numeric(survdata$OS),riskScore=xx)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(merge_data, time = "time", event = "event",
   variables = c("riskScore"))
summary(res.cut)
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
# 4. Fit survival curves and visualize
library("survival")
fit <- survfit(Surv(time, event) ~riskScore, data = res.cat)
gg<-ggsurvplot(fit,pval = TRUE, conf.int = TRUE,pval.coord = c(100, 0.6),
          ,ylim = c(0.1, 1), xlab = "Time(Days)",
		   risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
            ggtheme = theme_bw() # Change ggplot2 theme
           )
pdf(file="Fig.4A.pdf",height=8,width=10)
	print(gg, newpage = FALSE)
dev.off()

#####################METABRIC validation######################################
#####################METABRIC 验证######################################
#####################METABRIC 验证######################################
METABRIC_expr<-read.table("brca_metabric/data_expression_median.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
##########基因注释信息####	
METABRIC_LumA<-METABRIC_expr[METABRIC_expr[,1] %in% genes,] 
rownames(METABRIC_LumA)<-METABRIC_LumA[,1];METABRIC_LumA<-METABRIC_LumA[,-c(1,2)]
############临床数据########
clinData<-read.table("brca_metabric_clinical_data -subset.txt",sep="\t",header=TRUE,stringsAsFactor=FALSE);
clinData<-clinData[clinData[,3] %in% c("LumA"),]
clinData<-clinData[,c(1,4,6)]
head(clinData)
colnames(clinData)=c("samples","times","status")
# clinData[,1]<-chartr("-",".",clinData[,1]) 
clinData[clinData=="Living"]=0;clinData[clinData=="Died of Other Causes"]=0; ##censored data###
clinData[clinData=="Died of Disease"]=1
rownames(clinData)<-gsub("-", ".",clinData[,1]) ;
################  order sample###################
common_samples<-intersect(rownames(clinData),colnames(METABRIC_LumA))
a=0;b=0;
for(i in 1:length(common_samples)){ ##
	a[i]=which(rownames(clinData)==common_samples[i])  ##clinData
	b[i]=which(colnames(METABRIC_LumA)==common_samples[i])##profileData
}
clinData<-clinData[a,]
METABRIC_LumA<-METABRIC_LumA[,b]     ##
d=0;
for(i in 1:length(gene_order)){   ##
	d[i]=which(rownames(METABRIC_LumA)==gene_order[i])  ##
}
METABRIC_LumA<-METABRIC_LumA[d,]     ##
##########################################
xx<-c()
for (i in 1:ncol(METABRIC_LumA)){
	xx_temp<-0;
	for(j in 1:length(rr)){
		xx_temp<-xx_temp+rr[j]*as.numeric(METABRIC_LumA[j,i])
	}
	xx<-c(xx,xx_temp)
}
median<-median(as.numeric(xx))
clas<-ifelse(as.numeric(xx) > median ,"High-risk","Low-risk")
survdata<-cbind(clinData,class=clas)
###############
merge_data<-data.frame(time=as.numeric(survdata$times),event=as.numeric(survdata$status),riskScore=xx)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(merge_data, time = "time", event = "event",
   variables = c("riskScore"))
summary(res.cut)
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
# 4. Fit survival curves and visualize
library("survival")
fit <- survfit(Surv(time, event) ~riskScore, data = res.cat)
gg<-ggsurvplot(fit,pval = TRUE, conf.int = TRUE,pval.coord = c(100, 0.6),
          ,ylim = c(0.25, 1), xlab = "Time(Months)",
		   risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
            ggtheme = theme_bw() # Change ggplot2 theme
           )

pdf(file="Fig.4B.pdf",height=8,width=10)
	print(gg, newpage = FALSE)
dev.off()


############DEGs in Tumor and normal
############DEGs in Tumor and normal
genes<-c("E2F1")
###########expr mat###########
expr     <-read.table("TCGA-BRCA.htseq_fpkm.tsv",stringsAsFactors=FALSE,sep="\t")
colnames(expr)<-expr[1,];expr<-expr[-1,]
rownames(expr)<-expr[,1];expr<-expr[,-1]
anno_data<-read.table("gencode.v22.annotation.gene.probeMap",header=TRUE,sep="\t",stringsAsFactors=FALSE)
anno_data<-anno_data[anno_data$gene %in% genes,] 
########
for(i in 1:nrow(anno_data)){
	tmp_expr<-t(expr[which(rownames(expr)==anno_data[i,1]),])
	Tumor_sample<-grep(".01A$",rownames(tmp_expr),value=TRUE)
	Tumor_expr<-tmp_expr[rownames(tmp_expr) %in% Tumor_sample,]  ###select the tumor samples##
	Normal_sample<-grep(".11A$",rownames(tmp_expr),value=TRUE)
	Normal_expr<-tmp_expr[rownames(tmp_expr) %in% Normal_sample,]  ###select the tumor samples##
	plt_data<-data.frame(Group=c(rep("Tumor",length(Tumor_expr)),rep("Normal",length(Normal_expr))),
						tmp_expr=c( as.numeric(Tumor_expr), as.numeric(Normal_expr)))
	title_text<-paste("boxplot of ",anno_data$gene[i],sep="")

	gg<-ggplot(plt_data, aes(x=Group, y=log(tmp_expr+0.01))) + 
		ylab("Expression level (log2)")+ 
		geom_boxplot(
			color=c("#529D3F","#D36762"),
			fill =c("#A8CE9E","#E29C98"),		na.rm = TRUE
		)+
		labs(title=title_text)+
		stat_compare_means(label.x = 1.4)+
		theme_classic()
		ggsave(gg,file="Fig.4C.pdf",height=8,width=10)
}



BRCAsubtype<-read.table("BRCA_subtype.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
BRCAsubtype$AJCC.Stage[grepl("Stage IV", BRCAsubtype$AJCC.Stage)]<-"Stage 4";  ##we need to reverse
BRCAsubtype$AJCC.Stage[grepl("Stage III",BRCAsubtype$AJCC.Stage)]<-"Stage 3"; 	##or replaced byStage I
BRCAsubtype$AJCC.Stage[grepl("Stage II", BRCAsubtype$AJCC.Stage)]<-"Stage 2";
BRCAsubtype$AJCC.Stage[grepl("Stage I",  BRCAsubtype$AJCC.Stage)]<-"Stage 1";
BRCAsubtype$AJCC.Stage[grepl("Stage 1",  BRCAsubtype$AJCC.Stage)]<-"Stage I";
BRCAsubtype$AJCC.Stage[grepl("Stage 2",  BRCAsubtype$AJCC.Stage)]<-"Stage II";
BRCAsubtype$AJCC.Stage[grepl("Stage 3",  BRCAsubtype$AJCC.Stage)]<-"Stage III";
BRCAsubtype$AJCC.Stage[grepl("Stage 4",  BRCAsubtype$AJCC.Stage)]<-"Stage IV";
BRCAsubtype<-BRCAsubtype[BRCAsubtype$AJCC.Stage %in% c("Stage IV","Stage III","Stage II","Stage I"),]
rownames(BRCAsubtype)<-BRCAsubtype[,1]
###########expr mat############
expr     <-read.table("TCGA-BRCA.htseq_fpkm.tsv",stringsAsFactors=FALSE,sep="\t")
colnames(expr)<-expr[1,];expr<-expr[-1,]
rownames(expr)<-expr[,1];expr<-expr[,-1]
Tumor_sample<-grep(".01A$",colnames(expr),value=TRUE)
Tumor_expr<-expr[,colnames(expr) %in% Tumor_sample] 
anno_data<-read.table("gencode.v22.annotation.gene.probeMap",header=TRUE,sep="\t",stringsAsFactors=FALSE)
genes<-c("E2F1")
anno_data<-anno_data[anno_data$gene %in% genes,] 
for(i in 1:nrow(anno_data)){
	tmp_expr<-t(Tumor_expr[which(rownames(Tumor_expr)==anno_data[i,1]),])
	rowname<-gsub("-01A$", "",rownames(tmp_expr)) 
	rownames(tmp_expr)<-rowname;colnames(tmp_expr)<-"expr";
	merge_data <- merge(BRCAsubtype, tmp_expr, by=0) 

	plt_data<-data.frame(Stage=merge_data$AJCC.Stage,
						expr=as.numeric(merge_data$expr))
	title_text<-paste("violin plot of ",anno_data$gene[i],sep="")
	
	statistic<-kruskal.test(expr ~ Stage, data = plt_data) 
	statistic$p.value						
	gg <-ggplot(plt_data, aes(x=Stage, y=expr, fill=Stage)) + # fill=name allow to automatically dedicate a color for each group
	geom_violin(draw_quantiles = c( 0.5))+
	# geom_boxplot(width=0.1, color="black", alpha=0.2)+
	ylab("Expression level")+ 
	theme_classic()+
	theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)) +     
	scale_fill_manual(values=c("#CD877E","#7B82B3","#91BF92","#916EAD"))+
	stat_compare_means(label.x = length(unique(plt_data$Stage))/2 + 0.5)+
	labs(title=title_text)
	ggsave(gg,file="Fig.4D.pdf",height=8,width=8)	
}


######################DEGs in sensitive and resistant
drug_response<-read.table("TCGA_drug.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)
drug_response<-drug_response[drug_response$cancer %in% c("BRCA"),]
drug_response$pharmaceutical_therapy_drug_name<-toupper(drug_response$pharmaceutical_therapy_drug_name)
drug_response$uniformName<-toupper(drug_response$uniformName)


AI<-c("ANASTRAZOLUM","ANASTROZOLE","ARIMIDEX","ANASTROZOLUM","ANASTRAZOLE","","EXEMESTANE","AROMATASE ","AROMASIN","","","LETROZOLE","FEMARA","","","TAMOXIFEN","TAMOXIPHENE","TAMOXIFEN","NOLVADEX","NOVADEX")
drug<-intersect(unique(drug_response$uniformName),AI)
drug_response<-drug_response[drug_response$pharmaceutical_therapy_type %in% c("Hormone Therapy"),]  
table(drug_response$treatment_best_response) 

###########expr mat###########
expr     <-read.table("TCGA-BRCA.htseq_fpkm.tsv",stringsAsFactors=FALSE,sep="\t")
colnames(expr)<-expr[1,];expr<-expr[-1,]
rownames(expr)<-expr[,1];expr<-expr[,-1]
anno_data<-read.table("gencode.v22.annotation.gene.probeMap",header=TRUE,sep="\t",stringsAsFactors=FALSE)
genes<-c("E2F1")
anno_data<-anno_data[anno_data$gene %in% genes,] 

for(i in 1:nrow(anno_data)){
	tmp_expr<-t(expr[which(rownames(expr)==anno_data[i,1]),])
	colname<-gsub("-01A$", "",colnames(expr)) 
	Resistance_sample<-drug_response[drug_response$treatment_best_response %in% c("Clinical Progressive Disease","Stable Disease"),] 
	Resistance_group<-intersect(colname,Resistance_sample$bcr_patient_barcode)
	Resistance_expr<-tmp_expr[colname %in% Resistance_group,]  ###select the Resistance samples##
	
	Sensitive_sample<-drug_response[drug_response$treatment_best_response %in% c("Complete Response"),] 
	Sensitive_group<-intersect(colname,Sensitive_sample$bcr_patient_barcode)
	Sensitive_expr<-tmp_expr[colname %in% Sensitive_group,]  ###select the Resistance samples##
	
	plt_data<-data.frame(Group=c(rep("Sensitive",length(Sensitive_expr)),rep("Resistance",length(Resistance_expr))),
						tmp_expr=c( as.numeric(Sensitive_expr), as.numeric(Resistance_expr)))
	title_text<-paste("boxplot of ",anno_data$gene[i],sep="")

	gg<-ggplot(plt_data, aes(x=Group, y=log(tmp_expr+0.01))) + 
		ylab("Expression level (log2)")+ 
		geom_boxplot(
			color=c("#529D3F","#D36762"),
			fill =c("#A8CE9E","#E29C98"),		na.rm = TRUE
		)+
		labs(title=title_text)+
		stat_compare_means(label.x = 1.4)+
		theme_classic()
	ggsave(gg,file="Fig.4E.pdf",height=8,width=10)
}




	
	












