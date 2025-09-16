setwd("E:\\组学数据分析\\CyTOF\\浙大3个项目_原始数据\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师")
library(cytofkit2)
antigen<-readRDS("E:\\组学数据分析\\CyTOF\\CytoF_antigen.rds")

set.seed(100)
####cytofkit数据整合了很多读文件 降维 聚类步骤 但是一次只能读一个fcs文件 ###这个是核心函数 相当与完成所有步骤
#############################
# res <- cytofkit(fcsFiles = basal[1], 
#                 markers = antigen$parameter, 
#                 projectName = 'cytofkit_test',
#                 transformMethod = "cytofAsinh", 
#                 mergeMethod = "all",
#                 fixedNum = 500,                                    ## set at 500 for faster run 每个样本挑选多少个细胞参与分析
#                 dimReductionMethod = "tsne",
#                 clusterMethods = c("Rphenograph", "ClusterX"),    ## accept multiple methods
#                 visualizationMethods = c("tsne", "pca"),          ## accept multiple methods
#                 progressionMethod = "isomap",
#                 clusterSampleSize = 500,
#                 resultDir = dir2,
#                 saveResults = TRUE, 
#                 saveObject = TRUE)
# res
# class(res)
#############################

dir="E:\\组学数据分析\\CyTOF\\浙大3个项目_原始数据\\3个项目_原始数据\\5_gate---圈门后\\1_16samples_CD45+\\费老师"
file <- list.files(dir ,pattern='.fcs$', full=TRUE)
file
#245810 是ctrl
basal<-file[grep("C002|C004|C005|C008|C010",file)]
#parameters <- as.character(read.table(paraFile, header = TRUE)[,1])
parameters<-antigen$parameter

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
#lapply(basal_cytof,function(x){colnames(x)})

###################
data_transformed <- cytof_exprsExtract(fcsFile = file,comp = FALSE,transformMethod = "cytofAsinh")


###(2)直接函数内置 整合 每个样本堆积提取2万个细胞
combined_data_transformed <- cytof_exprsMerge(fcsFiles = basal, comp=FALSE,
                                              transformMethod = "arcsinh",
                                              mergeMethod = "ceil",
                                              fixedNum = 20000)
dim(combined_data_transformed)
colnames(combined_data_transformed)
table(duplicated(colnames(combined_data_transformed)))


####2 antigen 对应
####2 提取细胞子集 （我们是每个样本随机提取10000个细胞）
data_transformed<-combined_data_transformed
#data_transformed_1k <- data_transformed[1:100, ]

###3 执行聚类 run PhenoGraph #method :Rphenograph ClusterX DensVM FlowSOM
cluster_PhenoGraph <- cytof_cluster(xdata = data_transformed, method = "Rphenograph", Rphenograph_k = 65)#K值越小 聚类个数越多

## run ClusterX
data_transformed_1k_tsne <- cytof_dimReduction(data=data_transformed, method = "tsne")
cluster_ClusterX <- cytof_cluster(ydata = data_transformed,  method="ClusterX")

## run FlowSOM with cluster number 15 降维时候使用的K值 再IL1R1中是K=60 HX论文：K=15
cluster_FlowSOM <- cytof_cluster(xdata = data_transformed, method = "FlowSOM", FlowSOM_k = 60)

## combine data
data_1k_all <- cbind(data_transformed, data_transformed_1k_tsne, 
                     PhenoGraph = cluster_PhenoGraph,  
                     FlowSOM=cluster_FlowSOM)
data_1k_all <- as.data.frame(data_1k_all)
dim(data_1k_all)

####4可视化
## PhenoGraph plot on tsne
cytof_clusterPlot(data=data_1k_all, xlab="tsne_1", ylab="tsne_2", 
                  cluster="PhenoGraph", sampleLabel = FALSE,labelSize = 2)+theme(plot.title = element_blank())


## PhenoGraph cluster heatmap
PhenoGraph_cluster_median <- aggregate(. ~ PhenoGraph, data = data_1k_all, median)
cytof_heatmap(PhenoGraph_cluster_median[, 2:37], baseName = "PhenoGraph Cluster Median")
