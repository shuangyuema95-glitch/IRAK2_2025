#对count数据归一化
#BiocManager::install("DESeq2")
library("DESeq2") 
library(dplyr)

#读入数据
LPS<-read.delim("clipboard",header=T);dim(LPS);rownames(LPS)<-LPS$Geneid;LPS$Geneid<-NULL
basal<-read.delim("clipboard",header=T);dim(basal);rownames(basal)<-basal$Geneid;basal$Geneid<-NULL
R848<-read.delim("clipboard",header=T);dim(R848);rownames(R848)<-R848$Geneid;R848$Geneid<-NULL

#构建对Count的归一化函数
NorRNA<-function(data,treat){#data就是上一步数据 treat是什么处理LPS还是R848

lie=readline()
if(lie=="no"){
  data<-data
}else{
  data<-data[,match(lie,colnames(data))]
}

treat<-as.character(treat)
database=data[,c(1:ncol(data))];dim(data)
head(database)
database=database[complete.cases(database),];dim(database)

##构建factor,切忌名称中间不能添加如“_“ 等符号，否则会报错Error in DESeqDataSet(se, design = design, ignoreRank) : variables in design formula cannot contain NA: condition
cs1<-length(grep("MR|ZCH|MPQ",colnames(database)));cs1
cs2<-(ncol(database)-cs1)
##condition <- factor(c("h0","h0","h0","h6","h6","h6"))#这里我只是列出了两个处理的数据，有多少数据写多少
condition<-factor(c(rep(treat,cs1),rep("control",cs2)))
##DEseq2均一化 DESeqDataSetFromMatrix 用于创建Deseq2的对象

##count的多种标准化方法
dds <- DESeqDataSetFromMatrix(database, DataFrame(condition), design= ~ condition)
dds <- DESeq(dds) 
#sizeFactors(dds)#rlog标准化 

##归一化count数据 resdata是差异表达的结果列 后面是归一化的表达矩阵
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, 
                                                          normalized=TRUE)),by="row.names",sort=FALSE)
return(resdata)

# merge_list <- data.frame(database,resdata)
# head(merge_list)
# resdata <- merge_list
# head(resdata)

# ##2.差异数据保存
# resultsNames(dds)
# # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
# table(res$padj<0.05) #取P值小于0.05的结果
# res <- res[order(res$padj),]
# diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
# diff_gene_deseq2 <- row.names(diff_gene_deseq2)
# resdata1 <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
# 得到csv格式的差异表达分析结果
}
LPS_nor<-NorRNA(LPS,"LPS")
basal_nor<-NorRNA(basal,"patient")
R848_nor<-NorRNA(R848,'R848')


colnames(LPS_nor)
colnames(basal_nor)
colnames(R848_nor)

#暂存数据
setwd("C:\\Users\\35565\\Desktop\\组学数据分析\\IRAK2\\RNA_result")
save(LPS_nor,basal_nor,R848_nor,file="LPS_basal_NorCount.Rdata")

#绘制特定通路的热图
library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
ZC<-NFKB
NFKB<-NFKB2
heatPath<-function(data){
data1<-data[,-match(c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),colnames(data))]
nfkb<-data1[data1$Row.names%in%NFKB$gene,];table(duplicated(nfkb$Row.names))#没有重复
ifn1<-data1[data1$Row.names%in%IFN1$gene,];table(duplicated(ifn1$Row.names))#没有重复
ifn2<-data1[data1$Row.names%in%IFN2$gene,];table(duplicated(ifn2$Row.names))#没有重复
mapk<-data1[data1$Row.names%in%MAPK$gene,];table(duplicated(mapk$Row.names))#没有重复
rownames(nfkb)<-nfkb$Row.names;rownames(ifn1)<-ifn1$Row.names
rownames(ifn2)<-ifn2$Row.names;rownames(mapk)<-mapk$Row.names
nfkb$Row.names<-NULL;ifn1$Row.names<-NULL;ifn2$Row.names<-NULL;mapk$Row.names<-NULL
nfkb<-if(length(which(as.numeric(apply(nfkb,1,sum))==0))==0){nfkb}else{nfkb[-which(as.numeric(apply(nfkb,1,sum))==0),]}
ifn1<-if(length(which(as.numeric(apply(ifn1,1,sum))==0))==0){ifn1}else{ifn1[-which(as.numeric(apply(ifn1,1,sum))==0),]}
ifn2<-if(length(which(as.numeric(apply(ifn2,1,sum))==0))==0){ifn2}else{ifn2[-which(as.numeric(apply(ifn2,1,sum))==0),]}
mapk<-if(length(which(as.numeric(apply(mapk,1,sum))==0))==0){mapk}else{mapk[-which(as.numeric(apply(mapk,1,sum))==0),]}

#heatmap 函数需要去除某一行全部是0 的数据
bk = unique(c(seq(-2.5,2.5, length=100)))
nfkb_p<- pheatmap(nfkb,scale = "row",cluster_cols=F,cluster_rows = T,show_rownames = F,border_color=NA,,main = "IFN2", fontsize = 3,legend = T,treeheight_col = 0,treeheight_row = 0,
                  legend_breaks=seq(-2,2,0.5),color=colorRampPalette(c("Blue","White","Red"))(101),annotation_legend = T, breaks = bk,cellwidth = 8,cellheight = 1.8)                                     

ifn1_p<- pheatmap(ifn1,scale = "row",cluster_cols=F,cluster_rows = T,show_rownames = T,border_color=NA,,main = "IFN-I", fontsize = 6.5,legend = T,treeheight_col = 0,treeheight_row = 0,
                     legend_breaks=seq(-2,2,0.5),color=colorRampPalette(c("Blue","White","Red"))(101),annotation_legend = T, breaks = bk,cellwidth = 8,cellheight = 6.8)   
ifn2_p<- pheatmap(ifn2,scale = "row",cluster_cols=F,cluster_rows = T,show_rownames = T,border_color=NA,,main = "IFN-II", fontsize = 6.5,legend = T,treeheight_col = 0,treeheight_row = 0,
                  legend_breaks=seq(-2,2,0.5),color=colorRampPalette(c("Blue","White","Red"))(101),annotation_legend = T, breaks = bk,cellwidth = 8,cellheight = 6.8)   
mapk_p<-pheatmap(mapk,scale = "row",cluster_cols=F,cluster_rows = T,show_rownames = T,border_color=NA,,main = "MAPK", fontsize = 6.5,legend = T,treeheight_col = 0,treeheight_row = 0,
                 legend_breaks=seq(-2,2,0.5),color=colorRampPalette(c("Blue","White","Red"))(101),annotation_legend = T, breaks = bk,cellwidth = 8,cellheight = 6.8)   
}

#功能富集分析
library(fgsea)          # GSEA分析主程序
library(data.table)     # 数据处理
library(ggplot2)        # 画图处理
library(dplyr)          # 数据处理
library(msigdb)         # 包含基因集合，通常和GSEA分析共同使用
library(GSEABase)       # 可以提供GSEA基础结构和函数,也会被其他包调用
library(tidyverse)
library(dplyr)
library(Seurat)
library(openxlsx)
#####这一步比较考验网速
#msigdb.hs<-getMsigdb(org = 'hs',id = c("SYM", "EZID"))
#msigdb.hs  <-  appendKEGG(msigdb.hs)#追加KEGG集合 默认的MsigDB中没有KEGG数据集
#listCollections(msigdb.hs)
#提取MsigDB的HALLMARK
#save(msigdb.hs,file="gsea_MSIGDB.Rdata")
load("C:\\Users\\35565\\Desktop\\组学数据分析\\IRAK2\\RNA_result\\gsea_MSIGDB.Rdata")
load("C:\\Users\\35565\\Desktop\\组学数据分析\\IRAK2\\RNA_result\\1-NFKB热图结果.Rdata")
basal_nor
LPS_nor
R848_nor

library(tidyverse)
library(dplyr)
#####对list按照foldchange降序排序
hallmarks  <-  subsetCollection(msigdb.hs, 'h')
msigdb_ids <- geneIds(hallmarks)

LPS_nor[which(!is.na(LPS_nor$log2FoldChange)),]%>%dplyr::select(log2FoldChange,Row.names)%>%arrange(desc(log2FoldChange))->gs_LPS
LPS_list=as.numeric(gs_LPS$log2FoldChange);names(LPS_list)<-gs_LPS$Row.names

R848_nor[which(!is.na(R848_nor$log2FoldChange)),]%>%dplyr::select(log2FoldChange,Row.names)%>%arrange(desc(log2FoldChange))->gs_R848
R848_list=as.numeric(gs_R848$log2FoldChange);names(R848_list)<-gs_R848$Row.names


lps_gsea <- fgsea(pathways = msigdb_ids, 
                  stats = LPS_list,
                  eps = 0.0,
                  minSize = 15,
                  maxSize = 500)

R848_gsea <- fgsea(pathways = msigdb_ids, 
                   stats = R848_list,
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)
##绘制富集图
#筛选几个top通路
topLPSUp <- lps_gsea[ES > 0][head(order(pval), n = 10), pathway]
topLPSDown <- lps_gsea[ES < 0][head(order(pval), n = 20), pathway]


topR848Up<- R848_gsea[ES > 0][head(order(pval), n = 10), pathway]
topR848Down <- R848_gsea[ES < 0][head(order(pval), n = 20), pathway]

library(ggplot2)
##(1) 两个条件刺激之后都在Interferon-α Interferon-γ中上调
topLPSUp#MYC WNT IFN-α IFN-γ中上调
topLPSDown[c(3,11,13,18)]#补体 凋亡信号通路 TNF-NFKB 糖酵解信号通路中下调

######上调通路
LPS_upPathway<-list()
for(i in topLPSUp[3:4]){
  LPS_upPathway[[i]]<-plotEnrichment(msigdb_ids[[i]],
                                     LPS_list,ticksSize = ) + labs(title=i)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(plot.title = element_text(size=7.8,face="bold"))
  theme(
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm"),
    axis.text.x=element_text(color="black",family="Arail",size=6), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(color="black",family="Arail",size=6), #设置y轴刻度标签的字体属性
    axis.title.y=element_text(family="Arail",size=5.5,face="bold"), #设置y轴的标题的字体属性
    axis.title.x=element_text(family="Arail",size=5.5,face="bold"), #设置x轴的标题的字体属性
    axis.line = element_line(color = "black", size = 1),
  )
}
LPS_upPathway[[1]]
LPS_upPathway[[2]]
plot_grid(plotlist = LPS_upPathway,ncol=1,rel_heights = c(1,1),align="hv")

####下调通路
LPS_downPathway<-list()
for(i in topLPSDown[c(3,11,13,18)]){
  LPS_downPathway[[i]]<-plotEnrichment(msigdb_ids[[i]],
                                       LPS_list) + labs(title=i)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(plot.title = element_text(size=7.8,face="bold"))
  theme(
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm"),
    axis.text.x=element_text(color="black",family="Arail",size=6), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(color="black",family="Arail",size=6), #设置y轴刻度标签的字体属性
    axis.title.y=element_text(family="Arail",size=5.5,face="bold"), #设置y轴的标题的字体属性
    axis.title.x=element_text(family="Arail",size=5.5,face="bold"), #设置x轴的标题的字体属性
    axis.line = element_line(color = "black", size = 1),
  )
}
LPS_downPathway[[1]]
LPS_downPathway[[2]]
LPS_downPathway[[3]]
LPS_downPathway[[4]]
library(cowplot)
plot_grid()

  
##(2) 下调通路：
topR848Up#PI3K-AKT UPR
topR848Down[c(8,10,12)]#补体 血管生成 炎症应答



