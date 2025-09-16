# setwd("C:\\Users\\35565\\Desktop\\组学数据分析\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师")
# setwd("E:\\组学数据分析\\CyTOF\\浙大3个项目_原始数据\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师")
# library(cytofkit2)
# 
# ###############1 读入FCS文件
# rm(list=ls())
# require(cytofWorkflow)
# p1='C:\\Users\\35565\\Desktop\\组学数据分析\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师'
# p1='E:\\组学数据分析\\CyTOF\\浙大3个项目_原始数据\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师'
# fs1=list.files(p1,'*fcs')
# fs1
# samp<-read.flowSet(files=fs1,path=p1);length(samp)#之后的分析 都是基于这个CyTOF对象
# samp#samp是一个list格式的文件 有几个fcs文件 这个list就存储几个对象
# exprs(samp[[1]])[1:6,1:5]#取出每个对象的基因表达谱
# colnames(exprs(samp[[1]]))
# 
# ###########cytof内置数据集
# ########
# # 内置数据集
# # library(HDCytoData)
# # fs<-Bodenmiller_BCR_XL_flowSet()#这是R包内置的cytof数据 对网速要求比较高
# # fs
# # lapply(1:10,function(x){dim(exprs(fs[[x]]))})
# 
# #内置数据集的其他下载方式
# # ehub <- ExperimentHub() # create ExperimentHub instance 
# # ehub <- query(ehub, "HDCytoData") # find HDCytoData datasets 
# # md <- as.data.frame(mcols(ehub)) # retrieve metadata table
# ##########################################
# 
# 
# ###############2 构建SingleCellExperiment对象 
# #文件1:fcs文件对象
# lapply(1:10, function(x) dim( exprs(samp[[x]]) ))#每个样本的细胞数量差别比较大
# fs<-samp
# #文件2:metadata:文件名 样本ID 条件 患者ID
# md<-read.delim("clipboard",header = T)
# head(data.frame(md))
# 
# #文件3:抗体信息文件
# panel<-read.delim("clipboard",header=T)
# head(data.frame(panel))
# 
# ###需要准备抗体信息文件 第一列是fcs文件列名 第二列是抗体标签 第三列是抗体类型
# # url<-"http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"# 这个cytof的panel的抗体信息表格
# # panel<-"PBMC8_panel_v3.xlsx"#这个R代码下载不下来的时候 去官网上下载
# # download.file(file.path(url,panel),destfile=panel,mode="wb")
# # panel<-read_excel(panel)
# # head(data.frame(panel))
# # library(CATALYST)
# sce<-prepData(fs,panel,md,features=panel$fcs_colname)
# #Error in prepData(fs, panel, md, features = panel$fcs_colname) : 
# #不是所有的panel[[panel_cols$channel]] %in% colnames(fs)都是TRUE
# all(panel$fcs_colname %in% colnames(fs))#这里是FALSE的话就是会报错
# 
# 
# 
# #############################从列名不和准备的antigen报错 就尝试cytofkit2的流程了 大同小异 
# #########################cytofkit流程
# #####1 安装cytokit
# BiocManager::install("cytofkit")
# if(!require(devtools)){
#   install.packages("devtools") # If not already installed
# }
# #package ‘cytofkit2’ is not available for Bioconductor version '3.19'
# devtools::install_github("JinmiaoChenLab/cytofkit2", dependencies=TRUE)
# #报错1 不能找到编译包
# pkgbuild::check_build_tools(debug = TRUE)
# #需要编译包解决
# #Rtools不能通过常规install.packages()命令进行安装，需要通过installr包进行安装
# install.packages("installr")
# install.packages("stringr")    ###依赖包
# library(stringr)
# library(installr)
# install.Rtools()

library(cytofkit2)
library(CATALYST)
library(dplyr)
library(flowCore)
library(ggcyto) 
library(ggplot2)
library(mvtnorm)
library(openCyto) 
library(patchwork)
library(reshape2)
library(pander)# BiocManager::install("pander")
library(tidyverse)
library(dplyr)


#####2 预处理流程 https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/cytofkit/inst/doc/cytofkit_workflow.html#pre-processing
#cytofkit2 一个简单的示例
##########

set.seed(100)
dir <- system.file('extdata',package='cytofkit2')
file <- list.files(dir ,pattern='.fcs$', full=TRUE)
parameters <- list.files(dir, pattern='.txt$', full=TRUE)


####cytofkit数据整合了很多读文件 降维 聚类步骤 但是一次只能读一个fcs文件
###########coreFunction
###############
setwd("E:\\IRAK2\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师")
folder="E:\\IRAK2\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师\\LPS/"
antigen<-readRDS("E:\\组学数据分析\\CyTOF\\CytoF_antigen.rds")
dir2="E:\\IRAK2\\Cytof_data"
out="G:\\cellranger学习笔记\\3个项目_原始数据\\RESULT"
res <- cytofkit(fcsFiles = folder, #这里写的是文件的路径 不是文件叫啥就可以读取了
                markers = antigen$parameter, 
                projectName = 'cytofkit_test',
                transformMethod = "arcsinh", 
                mergeMethod = "ceil",
                fixedNum = 10000,                                    ## set at 500 for faster run
                dimReductionMethod = "tsne",
                clusterMethods = "Rphenograph",    ## accept multiple methods
                visualizationMethods = "tsne",          ## accept multiple methods
                progressionMethod = "NULL",
                #clusterSampleSize = 500, #每个簇统一的个数
                saveResults =  out
)
res
class(res)

########################################################
########################################################
#####################数据准备
setwd("C:\\Users\\35565\\Desktop\\组学数据分析\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师")
setwd("E:\\IRAK2\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师/")
dir <- system.file('extdata',package='cytofkit2')
dir2="C:\\Users\\35565\\Desktop\\组学数据分析\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师"
dir2="E:\\IRAK2\\cellranger学习笔记\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师/"

file <- list.files(dir ,pattern='.fcs$', full=TRUE)
file
file2 <- list.files(dir2 ,pattern='.fcs$', full=TRUE)
file2

basal<-file2[grep("C002|C004|C005|C008",file2)]
samp_basal<-c("LL_basal","WJ_basal","ZCH_basal","MR_basal")
LPS<-file2[grep("C001|C003|C006|C007|C009|C010",file2)]
samp_LPS<-c("WJ_LPS","LL_LPS","ZCH_LPS","MR_LPS","FYD_LPS","YYH_LPS")
antigen<-readRDS("C:\\Users\\35565\\Desktop\\组学数据分析\\CYTOF\\CytoF_antigen.rds")
antigen<-readRDS("E:\\组学数据分析\\CyTOF\\CytoF_antigen.rds")

####################
#paraFile <- list.files(dir, pattern='.txt$', full=TRUE)
#antigen<-read.delim("clipboard",header=T)
#parameters <- as.character(read.table(paraFile, header = TRUE)[,1])

####1 读入数据 每个样本随机选择10000个细胞 使用反正弦双曲变化 转化数据 提取表达矩阵 ###一次读入一个 然后merge
###(1)用list存储
##################
# data_transformed <- cytof_exprsExtract(fcsFile = file,
#                                        comp = FALSE,
#                                        transformMethod = "cytofAsinh")
# 
# basal_cytof<-list()
# for(i in 1:length(basal)){
# data_transformed <- cytof_exprsExtract(fcsFile = basal[1], 
#                                        comp = FALSE, 
#                                        transformMethod = "cytofAsinh")
# basal_cytof[[i]]<-data_transformed
# }
# lapply(basal_cytof,function(x){dim(x)})
#########################################
#lapply(basal,function(x){colnames(x)})
###################
# data_transformed <- cytof_exprsExtract(fcsFile = file,comp = FALSE,transformMethod = "cytofAsinh")
# LPS<-file2[grep("C001|C003|C006|C007|C009|C010",file2)]
# samp_LPS<-c("WJ_LPS","LL_LPS","ZCH_LPS","MR_LPS","FYD_LPS","YYH_LPS")

###(2)直接函数内置 整合 每个样本堆积提取2万个细胞
LPS_combine <- cytof_exprsMerge(fcsFiles = LPS, comp=FALSE,
                                transformMethod = "arcsinh",
                                mergeMethod = "ceil",
                                fixedNum = 18000)
basal_combine <- cytof_exprsMerge(fcsFiles = basal, comp=FALSE,
                                  transformMethod = "arcsinh",
                                  mergeMethod = "ceil",
                                  fixedNum = 18000)
dim(LPS_combine)
dim(basal_combine)
####添加样本信息
# LPS_combine<-as.data.frame(LPS_combine)
# LPS_combine%>%dplyr::mutate(sample=dplyr::case_when(
#   length(grep("_C001",rownames(LPS_combine)))!=0 ~ 'WJ',length(grep("_C003",rownames(LPS_combine)))!=0 ~ 'LL',
#   length(grep("_C006",rownames(LPS_combine)))!=0 ~ 'ZCH',length(grep("_C007",rownames(LPS_combine)))!=0 ~ 'MR',
#   length(grep("_C009",rownames(LPS_combine)))!=0 ~ 'FYD',length(grep("_C010",rownames(LPS_combine)))!=0 ~ 'YYH'))->LPS_combine
# table(LPS_combine$sample)
sampLPS<-data.frame(ID=rownames(LPS_combine))
LPS_combine<-as.data.frame(LPS_combine)
samplel<-apply(sampLPS,1,function(x){
  if(length(grep("_C001",x))!=0){
    "WJ"
  }else if(length(grep("_C003",x))!=0){
    "LL"
  }else if(length(grep("_C006",x))!=0){
    "ZCH"
  }else if(length(grep("_C007",x))!=0){
    "MR"
  }else if(length(grep("_C009",x))!=0){
    "FYD"
  }else{
    "YYH"
  }
})
L1<-paste(rownames(LPS_combine),samplel,sep="_")
rownames(LPS_combine)<-L1

#rownames(basal_combine)<-gsub("_YYH","",rownames(basal_combine))
rownames(basal_combine)
sampbasal<-data.frame(ID=rownames(basal_combine))
basal_combine<-as.data.frame(basal_combine)
sampleb<-apply(sampbasal,1,function(x){
  if(length(grep("_C002",x))!=0){
    "LL"
  }else if(length(grep("_C004",x))!=0){
    "WJ"
  }else if(length(grep("_C005",x))!=0){
    "ZCH"
  }else{
    "MR"
  }
})
L2<-paste(rownames(basal_combine),sampleb,sep="_")
rownames(basal_combine)<-L2
table(duplicated(colnames(LPS_combine)))#没有重复 下面使用match匹配

####2 antigen 对应
#antigen<-read.delim("clipboard",header=T)
# dim(antigen)
# #需要和FCS文件对一下列名
# ma<-c()
# for(i in 1:nrow(antigen)){
#   ma<-c(ma,colnames(LPS_combine)[grep(antigen$Label[i],colnames(LPS_combine))])
# }
# length(ma)
# antigen$Label[1:6]
# antigen$parameter<-ma
# saveRDS(antigen,file="CytoF_antigen.rds")

####2 提取细胞子集 （我们是每个样本随机提取10000个细胞）
LPS_trans<-LPS_combine
basal_trans<-basal_combine

###3 执行聚类 run PhenoGraph #method :Rphenograph ClusterX DensVM FlowSOM
PhenoGra_LPS <- cytof_cluster(xdata = LPS_trans, method = "Rphenograph",Rphenograph_k=105)#K值越小 聚类个数越多
## run ClusterX
tsne_LPS <- cytof_dimReduction(data=LPS_trans, method = "tsne")
ClusterX_LPS <- cytof_cluster(ydata = tsne_LPS,  method="ClusterX")
## run FlowSOM with cluster number 15 降维时候使用的K值 再IL1R1中是K=60 HX论文：K=15
FlowSOM_LPS <- cytof_cluster(xdata = LPS_trans, method = "FlowSOM", FlowSOM_k = 12)
## combine data
LPS_all <- cbind(LPS_trans, tsne_LPS, 
                 PhenoGraph = PhenoGra_LPS, ClusterX=ClusterX_LPS, 
                 FlowSOM=FlowSOM_LPS)
LPS_all <- as.data.frame(LPS_all)
dim(LPS_all)
colnames(LPS_all)

###3 执行聚类 run PhenoGraph #method :Rphenograph ClusterX DensVM FlowSOM
PhenoGra_basal <- cytof_cluster(xdata = basal_trans, method = "Rphenograph",Rphenograph_k=105)
## run ClusterX
tsne_basal <- cytof_dimReduction(data=basal_trans, method = "tsne")
ClusterX_basal <- cytof_cluster(ydata = tsne_basal,  method="ClusterX")
## run FlowSOM with cluster number 15 降维时候使用的K值 再IL1R1中是K=60 HX论文：K=15
FlowSOM_basal <- cytof_cluster(xdata = basal_trans, method = "FlowSOM", FlowSOM_k = 12)
## combine data
basal_all <- cbind(basal_trans, tsne_basal, 
                   PhenoGraph = PhenoGra_basal, ClusterX=ClusterX_basal, 
                   FlowSOM=FlowSOM_basal)
basal_all <- as.data.frame(basal_all)



####4可视化
library(ggplot2)
## PhenoGraph plot on tsne
LPS_tsne<-cytof_clusterPlot(data=LPS_all, xlab="tsne_1", ylab="tsne_2", 
                      cluster="PhenoGraph", sampleLabel = FALSE,addLabel = T,
                      point_size = 0.1)+theme(
    plot.title=element_blank(),
    legend.position="NONE", #不需要图例
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_text(size=10), #设置y轴的标题的字体属性
    axis.title.x=element_text(size=10), #设置x轴的标题的字体属性
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank()
  )
LPS_tsne 

## PhenoGraph cluster heatmap
LPS_PheGmedian <- aggregate(. ~PhenoGra_LPS, data = LPS_all, median)
cytof_heatmap(LPS_PheGmedian[,2:64] , baseName = "",scaleMethod = "none",
              dendrogram="none",cex_row_label = 0.3,cex_col_label = 0.24)

kangT<-match(antigen$parameter,colnames(LPS_PheGmedian))
Clu<-match("PhenoGra_LPS",colnames(LPS_PheGmedian))
LPS_marker5<-t(apply(LPS_PheGmedian,1,function(x){
  max5<-sort(as.numeric(x[kangT]),decreasing = TRUE)[1:10]
  return(c(as.numeric(x[1]),colnames(LPS_PheGmedian)[match(max5,x)]))
}))
colnames(LPS_marker5)<-c("cluster",paste(rep("antigen",10),1:10,sep="_"))
LPS_marker5<-as.data.frame(LPS_marker5)
LPS_marker5$antigen10<-apply(LPS_marker5,1,function(x){
  y<-as.character(x)[-1]
  chai<-c()
  for(i in 1:length((y))){
    chai<-c(chai,gsub(">","",unlist(strsplit(unlist(strsplit(y[i],split = "<"))[2],split="_"))[2]))
  }
  paste(chai,collapse=",")
})
library(openxlsx)
write.xlsx(LPS_marker5,file="LPS_10.xlsx")


library(openxlsx)
setwd("E:\\IRAK2\\Cytof_data")
#write.xlsx(LPS_marker5,file="LPS_marker.xlsx")



Marker<-readRDS("E:\\组学数据分析\\单细胞\\下游数据\\生信宝库_单细胞注释marker\\Marker.rds")
###########################################注释细胞群
#https://mp.weixin.qq.com/s/npLxjXoEWwCuppn6_cR-FQ
#https://www.sohu.com/a/294932205_610701
#https://www.sohu.com/a/437675430_610701
#https://www.cellsignal.cn/pathways/immune-cell-markers-human
#https://www.jianshu.com/p/0c745b965620
LPS_all%>%mutate(celltype=case_when(
  1 ~ '',
  
))

#########根据比例B 少 Treg少 CD8T 增加 确定一下这几个群




############################################
######basal状态的
basal_tsne<-cytof_clusterPlot(data=basal_all, xlab="tsne_1", ylab="tsne_2", 
                              cluster="PhenoGraph", sampleLabel = FALSE,addLabel = T,
                              point_size = 0.1)+theme(
    plot.title=element_blank(),
    legend.position="none", #不需要图例
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_text(size=10), #设置y轴的标题的字体属性
    axis.title.x=element_text(size=10), #设置x轴的标题的字体属性
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank()
  )
basal_tsne 

## PhenoGraph cluster heatmap
basal_PheGmedian <- aggregate(. ~PhenoGra_basal, data = basal_all, median)
cytof_heatmap(basal_PheGmedian , baseName = "PhenoGraph Cluster Median")


kangT<-match(antigen$parameter,colnames(basal_PheGmedian))
Clu<-match("PhenoGra_basal",colnames(basal_PheGmedian))
basal_marker5<-t(apply(basal_PheGmedian,1,function(x){
  max5<-sort(as.numeric(x[kangT]),decreasing = TRUE)[1:5]
  return(c(as.numeric(x[1]),colnames(basal_PheGmedian)[match(max5,x)]))
}))
colnames(basal_marker5)<-c("cluster",paste(rep("antigen",5),1:5,sep="_"))
class(basal_marker5)
library(openxlsx)
setwd("E:\\IRAK2\\Cytof_data")
write.xlsx(basal_marker5,file="basal_marker.xlsx")

##############################################################
##############################################################
###########Merge 
##################(1) GM-CSF在ctrl和LPS中都画一下图 数据不服从正太分布 用柱状图看一下
# GMc<- grep("GM_CSF",colnames(basal_trans))
# GMc
# GMc<- grep("GM_CSF",colnames(LPS_trans))
# GMc
# 
# GM_dat<-data.frame()
# samp<-c("LL","WJ","ZCH","MR")
# Heng<-LETTERS[1:8]
# for(i in 1:length(samp)){
# sha<-basal_trans[grep(samp[i],rownames(basal_trans)),'Er168Di<168Er_GM_CSF>']
# xia<-LPS_trans[grep(samp[i],rownames(LPS_trans)),'Er168Di<168Er_GM_CSF>']
# he<-c(rep(Heng[i],length(sha)),rep(Heng[2*i],length(xia)))
# condition<-c(rep("basal",length(sha)),rep("LPS",length(xia)))
# sample<-rep(samp[i],length(condition))
# GM_dat<-rbind(GM_dat,data.frame(sample=sample,he=he,exp=c(as.numeric(sha),as.numeric(xia)),condition=condition,heng=rep(i,length(condition))))
# }
# dim(GM_dat)

dim(GM_dat)
Gmp<-c()
for(i in 1:length(samp)){
  GM_dat%>%dplyr::filter(sample==samp[i])%>%dplyr::filter(condition=="basal")->y1
  GM_dat%>%dplyr::filter(sample==samp[i])%>%dplyr::filter(condition=="LPS")->y2
  Gmp<-c(Gmp,wilcox.test(as.numeric(y1$exp),as.numeric(y2$exp),paired = FALSE)$p.value)
}


########最大值
##########
se<-function(x){as.numeric(sd(x))/sqrt(length(x))}
samp<-c("LL","WJ","ZCH","MR")
max_GM<-data.frame()
for(i in 1:length(samp)){
sha<-as.numeric(basal_trans[grep(samp[i],rownames(basal_trans)),'Er168Di<168Er_GM_CSF>'])
xia<-as.numeric(LPS_trans[grep(samp[i],rownames(LPS_trans)),'Er168Di<168Er_GM_CSF>'])
s1=c(max(sha),max(xia))
s1
e1=c(se(sha),se(xia))
sample=rep(samp[i],length(s1))
condition=c("basal","LPS")
dat1<-data.frame(max=s1,se=e1,sample=sample,condition=condition)
max_GM<-rbind(max_GM,dat1)  
}
###########
max_GM$heng<-c("a","a","b","b","c","c","d","d")#####用数字的时候 添加不上去
max_GM_bar<-ggplot(max_GM, aes(heng,max, fill = condition)) +
  geom_col(position = position_dodge(width = 0.45), width = 0.45) +
  geom_errorbar(aes(ymin = max -20*se, ymax = max +20*se), width = 0.17, size = 0.8, position = position_dodge(0.45)) +
  scale_fill_manual(values =c("#204FA1","#E04227"))+
  labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
  scale_x_discrete(label=c("LL","WJ","ZCH","MR"))+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_classic()+theme(legend.position = 0)+
  theme_classic()+theme(legend.position = 0)+
  theme_bw()+ylab("expression (GM-CSF)")+
  #labs(color="NMF cluster")+
  theme(
    legend.key.height = unit(0.55, 'cm'),
    legend.key.width = unit(0.55, 'cm'),
    legend.position="top", #不需要图例
    axis.text.x=element_text(colour="black",size=10), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(colour="black",size=10), #设置y轴刻度标签的字体属性
    axis.title.y=element_text(size=12,face="bold"), #设置y轴的标题的字体属性
    axis.title.x=element_text(size=12,face="bold"), #设置x轴的标题的字体属性
    plot.title = element_text(size=12,face="bold",hjust = 0.5,lineheight=0.2), #设置总标题的字体属性
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank()
  )


############平均值
se<-function(x){as.numeric(sd(x))/sqrt(length(x))}
samp<-c("LL","WJ","ZCH","MR")
mean_GM<-data.frame()
for(i in 1:length(samp)){
  sha<-as.numeric(basal_trans[grep(samp[i],rownames(basal_trans)),'Er168Di<168Er_GM_CSF>'])
  xia<-as.numeric(LPS_trans[grep(samp[i],rownames(LPS_trans)),'Er168Di<168Er_GM_CSF>'])
  s1=c(mean(sha),mean(xia))
  s1
  e1=c(se(sha),se(xia))
  sample=rep(samp[i],length(s1))
  condition=c("basal","LPS")
  dat1<-data.frame(mean=s1,se=e1,sample=sample,condition=condition)
  mean_GM<-rbind(mean_GM,dat1)  
}
###########
mean_GM$heng<-c(1,1,2,2,3,3,4,4)
mean_GM_bar<-ggplot(mean_GM, aes(heng,mean, fill = condition)) +
  geom_col(position = position_dodge(width = 0.45), width = 0.45) +
  geom_errorbar(aes(ymin = mean -20*se, ymax = mean +20*se), width = 0.17, size = 0.8, position = position_dodge(0.45)) +
  scale_fill_manual(values =c("#204FA1","#E04227"))+
  labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_classic()+theme(legend.position = 0)+
  theme_classic()+theme(legend.position = 0)+
  theme_bw()+ylab("expression (GM-CSF)")+
  #labs(color="NMF cluster")+
  theme(
    legend.key.height = unit(0.55, 'cm'),
    legend.key.width = unit(0.55, 'cm'),
    legend.position="top", #不需要图例
    axis.text.x=element_text(colour="black",size=10), #设置x轴刻度标签的字体属性
    axis.text.y=element_text(colour="black",size=10), #设置y轴刻度标签的字体属性
    axis.title.y=element_text(size=12,face="bold"), #设置y轴的标题的字体属性
    axis.title.x=element_text(size=12,face="bold"), #设置x轴的标题的字体属性
    plot.title = element_text(size=12,face="bold",hjust = 0.5,lineheight=0.2), #设置总标题的字体属性
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(), #不显示网格线
    panel.grid.minor = element_blank()
  )






