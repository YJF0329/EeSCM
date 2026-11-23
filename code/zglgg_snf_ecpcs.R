library(MASS)
library(SNFtool)
library("Matrix")
library("parallel")
library(survival)
library(survminer)
library("R.matlab")

#设置为当前代码路径
setwd("F:\\课题代码")


#临床数据行是样本
KIRC_clinical=read.csv("F:\\课题代码\\MATLAB TCGA DATA\\CHONGZHENG\\中国LGG数据\\clinic_LGG_10ys.csv",header = TRUE)

####################################################################
#   method 6  SNF directly                                         #
#
#################################################################### 
#载入matlab数据
group2_snf = readMat("group2_snf.mat")
group2_snf=as.integer(group2_snf[["group2.snf"]])
group3_snf = readMat("group3_snf.mat")
group3_snf=as.integer(group3_snf[["group3.snf"]])
group4_snf = readMat("group4_snf.mat")
group4_snf=as.integer(group4_snf[["group4.snf"]])
group5_snf = readMat("group5_snf.mat")
group5_snf=as.integer(group5_snf[["group5.snf"]])
group6_snf = readMat("group6_snf.mat")
group6_snf=as.integer(group6_snf[["group6.snf"]])
group7_snf = readMat("group7_snf.mat")
group7_snf=as.integer(group7_snf[["group7.snf"]])
group8_snf = readMat("group8_snf.mat")
group8_snf=as.integer(group8_snf[["group8.snf"]])
group9_snf = readMat("group9_snf.mat")
group9_snf=as.integer(group9_snf[["group9.snf"]])
group10_snf = readMat("group10_snf.mat")
group10_snf=as.integer(group10_snf[["group10.snf"]])

KIRC_SNF_lable<-list()
KIRC_SNF_lable[[1]]=group2_snf
KIRC_SNF_lable[[2]]=group3_snf
KIRC_SNF_lable[[3]]=group4_snf
KIRC_SNF_lable[[4]]=group5_snf
KIRC_SNF_lable[[5]]=group6_snf
KIRC_SNF_lable[[6]]=group7_snf
KIRC_SNF_lable[[7]]=group8_snf
KIRC_SNF_lable[[8]]=group9_snf
KIRC_SNF_lable[[9]]=group10_snf
#####log-rank检验
p_snf =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_SNF_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_SNF_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_snf[[temp]]= p_value
}      


####################################################################
#
#   ##NE                                         #
#
#################################################################### 
group2_NE = readMat("group2_NE.mat")
group2_NE =as.integer(group2_NE[["group2.NE"]])
group3_NE = readMat("group3_NE.mat")
group3_NE =as.integer(group3_NE[["group3.NE"]])
group4_NE = readMat("group4_NE.mat")
group4_NE =as.integer(group4_NE[["group4.NE"]])
group5_NE = readMat("group5_NE.mat")
group5_NE =as.integer(group5_NE[["group5.NE"]])
group6_NE = readMat("group6_NE.mat")
group6_NE =as.integer(group6_NE[["group6.NE"]])
group7_NE = readMat("group7_NE.mat")
group7_NE =as.integer(group7_NE[["group7.NE"]])
group8_NE = readMat("group8_NE.mat")
group8_NE =as.integer(group8_NE[["group8.NE"]])
group9_NE = readMat("group9_NE.mat")
group9_NE =as.integer(group9_NE[["group9.NE"]])
group10_NE = readMat("group10_NE.mat")
group10_NE =as.integer(group10_NE[["group10.NE"]])
#2—10次聚类结果（原命令）KIRP_W
KIRC_NE_lable<-list()
KIRC_NE_lable[[1]]=group2_NE
KIRC_NE_lable[[2]]=group3_NE
KIRC_NE_lable[[3]]=group4_NE
KIRC_NE_lable[[4]]=group5_NE
KIRC_NE_lable[[5]]=group6_NE
KIRC_NE_lable[[6]]=group7_NE
KIRC_NE_lable[[7]]=group8_NE
KIRC_NE_lable[[8]]=group9_NE
KIRC_NE_lable[[9]]=group10_NE

#####log-rank检验
p_NE =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_NE_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_NE_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_NE[[temp]]= p_value
}  

####################################################################
#
#   ##skf                                         #
#
####################################################################

#2—10次聚类结果（原命令）KIRP_W

group2_skf = readMat("group2_skf.mat")
group2_skf =as.integer(group2_skf[["group2.skf"]])
group3_skf = readMat("group3_skf.mat")
group3_skf =as.integer(group3_skf[["group3.skf"]])
group4_skf = readMat("group4_skf.mat")
group4_skf =as.integer(group4_skf[["group4.skf"]])
group5_skf = readMat("group5_skf.mat")
group5_skf =as.integer(group5_skf[["group5.skf"]])
group6_skf = readMat("group6_skf.mat")
group6_skf =as.integer(group6_skf[["group6.skf"]])
group7_skf = readMat("group7_skf.mat")
group7_skf =as.integer(group7_skf[["group7.skf"]])
group8_skf = readMat("group8_skf.mat")
group8_skf =as.integer(group8_skf[["group8.skf"]])
group9_skf = readMat("group9_skf.mat")
group9_skf =as.integer(group9_skf[["group9.skf"]])
group10_skf = readMat("group10_skf.mat")
group10_skf =as.integer(group10_skf[["group10.skf"]])
KIRC_skf_lable<-list()
KIRC_skf_lable[[1]]=group2_skf
KIRC_skf_lable[[2]]=group3_skf
KIRC_skf_lable[[3]]=group4_skf
KIRC_skf_lable[[4]]=group5_skf
KIRC_skf_lable[[5]]=group6_skf
KIRC_skf_lable[[6]]=group7_skf
KIRC_skf_lable[[7]]=group8_skf
KIRC_skf_lable[[8]]=group9_skf
KIRC_skf_lable[[9]]=group10_skf
#####log-rank检验
p_skf =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_skf_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_skf_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_skf[[temp]]= p_value
} 



####################################################################
#
#   ##hopes                                         #
#
####################################################################
group2_hopes = readMat("group2_hopes.mat")
group2_hopes =as.integer(group2_hopes[["group2.hopes"]])
group3_hopes = readMat("group3_hopes.mat")
group3_hopes =as.integer(group3_hopes[["group3.hopes"]])
group4_hopes = readMat("group4_hopes.mat")
group4_hopes =as.integer(group4_hopes[["group4.hopes"]])
group5_hopes = readMat("group5_hopes.mat")
group5_hopes =as.integer(group5_hopes[["group5.hopes"]])
group6_hopes = readMat("group6_hopes.mat")
group6_hopes =as.integer(group6_hopes[["group6.hopes"]])
group7_hopes = readMat("group7_hopes.mat")
group7_hopes =as.integer(group7_hopes[["group7.hopes"]])
group8_hopes = readMat("group8_hopes.mat")
group8_hopes =as.integer(group8_hopes[["group8.hopes"]])
group9_hopes = readMat("group9_hopes.mat")
group9_hopes =as.integer(group9_hopes[["group9.hopes"]])
group10_hopes = readMat("group10_hopes.mat")
group10_hopes =as.integer(group10_hopes[["group10.hopes"]])

#2—10次聚类结果（原命令）KIRP_W
KIRC_hopes_lable<-list()
KIRC_hopes_lable[[1]]=group2_hopes
KIRC_hopes_lable[[2]]=group3_hopes
KIRC_hopes_lable[[3]]=group4_hopes
KIRC_hopes_lable[[4]]=group5_hopes
KIRC_hopes_lable[[5]]=group6_hopes
KIRC_hopes_lable[[6]]=group7_hopes
KIRC_hopes_lable[[7]]=group8_hopes
KIRC_hopes_lable[[8]]=group9_hopes
KIRC_hopes_lable[[9]]=group10_hopes
#####log-rank检验
p_hopes =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_hopes_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_hopes_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_hopes[[temp]]= p_value
} 

####################################################################
#
#   ##MOVICS                
#
#
####################################################################
W_snf = readMat("W_snf.mat")#列表
W_snf=W_snf[["W.snf"]]#矩阵
W_NE = readMat("W_NE.mat")
W_NE =W_NE[["W.NE"]]
W_skf = readMat("W_skf.mat")
W_skf =W_skf[["W.skf"]]
W_hopes = readMat("W_hopes.mat")
W_hopes =W_hopes[["W.hopes"]]


clusterk_snf = estimateNumberOfClustersGivenGraph(W_snf, 2:5)
clusterk_ne = estimateNumberOfClustersGivenGraph(W_NE, 2:5)
clusterk_skf = estimateNumberOfClustersGivenGraph(W_skf, 2:5)
clusterk_hopes = estimateNumberOfClustersGivenGraph(W_hopes, 2:5)
#clusterk_movics = estimateNumberOfClustersGivenGraph(similarity.matrix, 2:5)
#MOVICS聚类数为2

#clusterk_snf = estimateNumberOfClustersGivenGraph(W_snf, 2:10)
#clusterk_ne = estimateNumberOfClustersGivenGraph(W_NE, 2:10)
#clusterk_skf = estimateNumberOfClustersGivenGraph(W_skf, 2:10)
#clusterk_hopes = estimateNumberOfClustersGivenGraph(W_hopes, 2:10)

#MOVICS聚类数为2
library(MOVICS)
moic1.res.list=NULL
samID=row.names(KIRC_clinical)
moic1.res.list[["SNF"]][["fit"]]=W_snf
clust.res.snf=data.frame(samID,clust=group5_snf)
moic1.res.list[["SNF"]][["clust.res"]]=clust.res.snf
moic1.res.list[["SNF"]][["mo.method"]]="SNF"

moic1.res.list[["NE"]][["fit"]]=W_NE
clust.res.NE=data.frame(samID,clust=group4_NE)
moic1.res.list[["NE"]][["clust.res"]]=clust.res.NE
moic1.res.list[["NE"]][["mo.method"]]="NE"

moic1.res.list[["skf"]][["fit"]]=W_skf
clust.res.skf=data.frame(samID,clust=group4_skf)
moic1.res.list[["skf"]][["clust.res"]]=clust.res.skf
moic1.res.list[["skf"]][["mo.method"]]="skf"

moic1.res.list[["hopes"]][["fit"]]=W_hopes
clust.res.hopes=data.frame(samID,clust=group2_hopes)
moic1.res.list[["hopes"]][["clust.res"]]=clust.res.hopes
moic1.res.list[["hopes"]][["mo.method"]]="hopes"

moic.res.list=moic1.res.list
moic.res <- lapply(moic.res.list, function(x) x$clust.res)
N.clust <- 2
fixed.sam <- moic.res[[1]]$samID

consensus.matrix <- matrix(0, nrow=length(fixed.sam), ncol=length(fixed.sam))
for (i in 1:length(moic.res)) {
  res <- moic.res[[i]]
  grpInfo <- res$clust; names(grpInfo) <- res$samID
  
  ans = as.data.frame(matrix(0, nrow=length(fixed.sam), ncol=length(fixed.sam)))
  rownames(ans) = colnames(ans) <- fixed.sam
  unique.grp = unique(grpInfo)
  
  for (k in 1:length(unique.grp)){
    grpK.members = names(grpInfo)[grpInfo==unique.grp[k]]
    ans[grpK.members, grpK.members] = 1
  }
  consensus.matrix <- consensus.matrix + as.matrix(ans)
}

similarity.matrix <- as.data.frame(consensus.matrix/length(moic.res))
# hcs <- hclust(vegdist(as.matrix(1-similarity.matrix), method = "jaccard"), "ward.D")
distance      = "euclidean"
#原方法，将linkage统一为average
hcs <- hclust(ClassDiscovery::distanceMatrix(as.matrix(1-similarity.matrix), distance), "average")
group2_movics <- cutree(hcs, N.clust)#集成后的聚类
group2_movics=as.integer(group2_movics)



#MOVICS聚类数为3

N.clust <- 3

group3_movics <- cutree(hcs, N.clust)#集成后的聚类
group3_movics=as.integer(group3_movics)


#MOVICS聚类数为4

N.clust <- 4

group4_movics <- cutree(hcs, N.clust)#集成后的聚类
group4_movics=as.integer(group4_movics)


#MOVICS聚类数为5

N.clust <- 5

group5_movics <- cutree(hcs, N.clust)#集成后的聚类
group5_movics=as.integer(group5_movics)


#MOVICS聚类数为6

N.clust <- 6

group6_movics <- cutree(hcs, N.clust)#集成后的聚类
group6_movics=as.integer(group6_movics)

#MOVICS聚类数为7

N.clust <- 7

group7_movics <- cutree(hcs, N.clust)#集成后的聚类
group7_movics=as.integer(group7_movics)


#MOVICS聚类数为8

N.clust <- 8

group8_movics <- cutree(hcs, N.clust)#集成后的聚类
group8_movics=as.integer(group8_movics)


#MOVICS聚类数为9

N.clust <- 9

group9_movics <- cutree(hcs, N.clust)#集成后的聚类
group9_movics=as.integer(group9_movics)


#MOVICS聚类数为10

N.clust <- 10

group10_movics <- cutree(hcs, N.clust)#集成后的聚类
group10_movics=as.integer(group10_movics)
#2—10次聚类结果（原命令）KIRP_W
KIRC_movics_lable<-list()
KIRC_movics_lable[[1]]=group2_movics
KIRC_movics_lable[[2]]=group3_movics
KIRC_movics_lable[[3]]=group4_movics
KIRC_movics_lable[[4]]=group5_movics
KIRC_movics_lable[[5]]=group6_movics
KIRC_movics_lable[[6]]=group7_movics
KIRC_movics_lable[[7]]=group8_movics
KIRC_movics_lable[[8]]=group9_movics
KIRC_movics_lable[[9]]=group10_movics
#####log-rank检验
p_movics =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_movics_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_movics_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_movics[[temp]]= p_value
} 

clusterk_movics = estimateNumberOfClustersGivenGraph(as.matrix(similarity.matrix), 2:5)
clusterk_movics2_10 = estimateNumberOfClustersGivenGraph(as.matrix(similarity.matrix), 2:10)
ECA = readMat("ECA.mat")#列表
ECA=ECA[["ECA"]]#矩阵
clusterk_ECA = estimateNumberOfClustersGivenGraph(ECA, 2:5)
clusterk_ECA2_10 = estimateNumberOfClustersGivenGraph(ECA, 2:10)
####################################################################
#
#   ##ECPCS                                       #
#
####################################################################
group2_ECPCS_HC = readMat("group2_ECPCS_HC.mat")
group2_ECPCS_HC =as.integer(group2_ECPCS_HC[["group2.ECPCS.HC"]])
group3_ECPCS_HC = readMat("group3_ECPCS_HC.mat")
group3_ECPCS_HC =as.integer(group3_ECPCS_HC[["group3.ECPCS.HC"]])
group4_ECPCS_HC = readMat("group4_ECPCS_HC.mat")
group4_ECPCS_HC =as.integer(group4_ECPCS_HC[["group4.ECPCS.HC"]])
group5_ECPCS_HC = readMat("group5_ECPCS_HC.mat")
group5_ECPCS_HC =as.integer(group5_ECPCS_HC[["group5.ECPCS.HC"]])
group6_ECPCS_HC = readMat("group6_ECPCS_HC.mat")
group6_ECPCS_HC =as.integer(group6_ECPCS_HC[["group6.ECPCS.HC"]])
group7_ECPCS_HC = readMat("group7_ECPCS_HC.mat")
group7_ECPCS_HC =as.integer(group7_ECPCS_HC[["group7.ECPCS.HC"]])
group8_ECPCS_HC = readMat("group8_ECPCS_HC.mat")
group8_ECPCS_HC =as.integer(group8_ECPCS_HC[["group8.ECPCS.HC"]])
group9_ECPCS_HC = readMat("group9_ECPCS_HC.mat")
group9_ECPCS_HC =as.integer(group9_ECPCS_HC[["group9.ECPCS.HC"]])
group10_ECPCS_HC = readMat("group10_ECPCS_HC.mat")
group10_ECPCS_HC =as.integer(group10_ECPCS_HC[["group10.ECPCS.HC"]])
#2—10次聚类结果（原命令）KIRP_W
KIRC_ECPCS_lable<-list()
KIRC_ECPCS_lable[[1]]=group2_ECPCS_HC
KIRC_ECPCS_lable[[2]]=group3_ECPCS_HC
KIRC_ECPCS_lable[[3]]=group4_ECPCS_HC
KIRC_ECPCS_lable[[4]]=group5_ECPCS_HC
KIRC_ECPCS_lable[[5]]=group6_ECPCS_HC
KIRC_ECPCS_lable[[6]]=group7_ECPCS_HC
KIRC_ECPCS_lable[[7]]=group8_ECPCS_HC
KIRC_ECPCS_lable[[8]]=group9_ECPCS_HC
KIRC_ECPCS_lable[[9]]=group10_ECPCS_HC

#####log-rank检验
p_ECPCS =matrix(NA,9)
temp=0
for(i in 1:9){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_ECPCS_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_ECPCS_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_ECPCS[[temp]]= p_value
} 


P=matrix(NA,nrow=6,ncol=9)
P[1,]=p_snf
P[2,]=p_NE
P[3,]=p_skf
P[4,]=p_hopes
P[5,]=p_movics
P[6,]=p_ECPCS
write.csv(P,"P.csv")







setwd("F:\\ecpcs结果\\LGG下游")
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ KIRC_ECPCS_lable[[1]], 
             data = KIRC_surv)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_value <- signif(p_value,3)
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
tiff(filename = "LGG.tiff",res = 300,width = 5,height = 3.8,units = 'in',compression = c("lzw"))
par(family="myFont") #设置后图里面字体均为Times New Roman
fit<-survfit(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_ECPCS_lable[[1]],data=KIRC_surv)
ggsurvplot(fit, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#FF5F5F", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_value, # 显示 pvalue
           pval.size = 6
)

dev.off()


group_sparse_kmeans=KIRC_ECPCS_lable[[1]]
#################################################################################################
#
# 寻找差异基因
#
################################################################################################


group_sample_sort=sort(group_sparse_kmeans,decreasing = FALSE,index.return=T)
lable_groupsample=group_sample_sort$ix
group_sample=data.frame(as.character(group_sample_sort$x))# sample分组文件
rownames(KIRC_clinical)=KIRC_clinical[,2]
KIRC_clinical=KIRC_clinical[,-c(1,2)]
rownames(group_sample)=rownames(KIRC_clinical)[lable_groupsample] #将排序后的分组文件行名也按照分组排列顺序重新排列
colnames(group_sample)='Cluster'
#write.csv(group_sample,"group_sample.csv")

#colors=list(Cluster=c('1'="#ff7473", '2'="#4bc2c5",'3'="#ffc33c"))
colors=list(Cluster=c('1'="#ff7473",'2'="#4bc2c5"))

###############################################################################
# location for each cluster #
###############################################################################
table(group_sparse_kmeans)

location_cluster1=which(group_sparse_kmeans==1)
location_cluster2=which(group_sparse_kmeans==2)
#location_cluster3=which(group_sparse_kmeans==3)
#location_cluster4=which(group_sparse_kmeans==4)

group_sparse_kmeans <- as.matrix(group_sparse_kmeans)
group_sparse_kmeans <- as.integer(group_sparse_kmeans)
## 定义下调上调表达量
count_under=function(x){
  length(which(x==-1))
}

count_over=function(x){
  length(which(x==1))
}

miRNA_LGG_delmeancv=read.csv("F:\\课题代码\\MATLAB TCGA DATA\\CHONGZHENG\\中国LGG数据\\quant10 3omics\\miRNA_LGG_delmeancv.csv",header=T,stringsAsFactors = F,row.names = 1)
mRNA_LGG_delmeancv=read.csv("F:\\课题代码\\MATLAB TCGA DATA\\CHONGZHENG\\中国LGG数据\\quant10 3omics\\mRNA_LGG_delmeancv.csv",header=T,stringsAsFactors = F,row.names = 1)
methy_LGG_delmeancv=read.csv("F:\\课题代码\\MATLAB TCGA DATA\\CHONGZHENG\\中国LGG数据\\quant10 3omics\\methy_LGG_delmeancv.csv",header=T,stringsAsFactors = F,row.names = 1)

SmiRNA_KIRP=standardNormalization(miRNA_LGG_delmeancv)
SmRNA_KIRP=standardNormalization(mRNA_LGG_delmeancv)
Smethy_KIRP=standardNormalization(methy_LGG_delmeancv)
########################################################################################################
# 4.1 miRNA differential analysis  #
#
#########################################################################################################

#quantile(SmiRNA_KIRP, probs = c(0.25,0.75))
#quantile(SmiRNA_KIRP)
setwd("F:\\ecpcs结果\\LGG下游\\differential gene/miRNA")

#######################################################################################
# step 1 kruskal.test for each miRNA   #
######################################################################################
p_diffmiRNA_all=NULL
for(i in 1:dim(miRNA_LGG_delmeancv)[2]){
  fit_kwH_miRNA_all=NULL
  fit_kwH_miRNA_all=kruskal.test(miRNA_LGG_delmeancv[,i]~as.factor(group_sparse_kmeans))
  p_diffmiRNA_all=c(p_diffmiRNA_all,fit_kwH_miRNA_all$p.value)
}

#######FDR adjust########
#using the p.adjust function
q_adjust_miRNA_all <- p.adjust (p_diffmiRNA_all, method = "fdr")
length(which(q_adjust_miRNA_all<0.05))
differ_miRNA_fdr=colnames(miRNA_LGG_delmeancv)[which(q_adjust_miRNA_all<0.05)]
differ_miRNA_fdr=as.data.frame(differ_miRNA_fdr)
write.csv(differ_miRNA_fdr,file='differ_miRNA_fdr.csv')
save(differ_miRNA_fdr,file='differ_miRNA_fdr.Rdata')



########################################################################################
#step 2 hypergeometric test for the fdar significant miRNA   #
########################################################################################
length(which(q_adjust_miRNA_all>=0.05))

#based on the scaled data,and choose the significant miRNA 定义<0.5为-1，>=0.5为1
bin_SmiRNA_KIRP=SmiRNA_KIRP[,which(q_adjust_miRNA_all<0.05)]
for(i in 1:dim(bin_SmiRNA_KIRP)[2]){
  for(j in 1:dim(bin_SmiRNA_KIRP)[1]){
    if(bin_SmiRNA_KIRP[j,i]<=-0.5){
      bin_SmiRNA_KIRP[j,i]=-1
    }else if(bin_SmiRNA_KIRP[j,i]>=0.5){
      bin_SmiRNA_KIRP[j,i]=1
    }else{
      bin_SmiRNA_KIRP[j,i]=0
    }
  }
}
save(bin_SmiRNA_KIRP,file='bin_SmiRNA_KIRP.Rdata')

hyper_miRNA_underover_final=list()
for (m in 1:2) {
  # under expressed genes z<=-1 for cluster
  location_cluster=get(paste('location_cluster',m,sep=''))
  #location_cluster=location_cluster2
  p_hyper_cluster_under=NULL
  for(i in 1:dim(bin_SmiRNA_KIRP)[2]){
    
    k_cluster=length(which(bin_SmiRNA_KIRP[location_cluster,i]==-1))
    m_cluster=length(which(bin_SmiRNA_KIRP[,i]==-1))
    N_cluster=dim(bin_SmiRNA_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_under=c(p_hyper_cluster_under,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_under_fdr =p.adjust(p_hyper_cluster_under, method = "fdr")
  p_hyper_cluster_under_fdr_sort=sort(p_hyper_cluster_under_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_under_fdr_sortp=p_hyper_cluster_under_fdr_sort$x
  hyper_fdr_cluster_miRNA_under=colnames(bin_SmiRNA_KIRP)[p_hyper_cluster_under_fdr_sort$ix][which(p_hyper_cluster_under_fdr_sortp<0.01)]
  hyper_fdr_cluster_miRNA_under=as.data.frame(hyper_fdr_cluster_miRNA_under)
  
  #for over expressed genes z>=1
  ##############################################
  # over expressed genes z>=1 for cluster  #
  ##############################################
  p_hyper_cluster_over=NULL
  for(i in 1:dim(bin_SmiRNA_KIRP)[2]){
    k_cluster=length(which(bin_SmiRNA_KIRP[location_cluster,i]==1))
    m_cluster=length(which(bin_SmiRNA_KIRP[,i]==1))
    N_cluster=dim(bin_SmiRNA_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_over=c(p_hyper_cluster_over,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_over_fdr= p.adjust(p_hyper_cluster_over, method = "fdr")
  p_hyper_cluster_over_fdr_sort=sort(p_hyper_cluster_over_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_over_fdr_sortp=p_hyper_cluster_over_fdr_sort$x
  hyper_fdr_cluster_miRNA_over=colnames(bin_SmiRNA_KIRP)[p_hyper_cluster_over_fdr_sort$ix][which(p_hyper_cluster_over_fdr_sortp <0.01)]
  hyper_fdr_cluster_miRNA_over=as.data.frame(hyper_fdr_cluster_miRNA_over)
  
  #combine under and over
  hyper_miRNA_underover=union(hyper_fdr_cluster_miRNA_under,
                              hyper_fdr_cluster_miRNA_over)
  hyper_miRNA_underover_final=c(hyper_miRNA_underover_final,(hyper_miRNA_underover))
  
}
save(hyper_miRNA_underover_final,file='hyper_miRNA_underover_final.Rdata')

#############combine all the under and over expression genes #############
hyper_miRNA_underover_final23=Reduce(union,hyper_miRNA_underover_final)
write.csv(hyper_miRNA_underover_final23,file='hyper_miRNA_underover_final23.csv',row.names = F)
save(hyper_miRNA_underover_final23,file='hyper_miRNA_underover_final23.Rdata')
#temp=match(hyper_miRNA_underover_final23,colnames(bin_SmiRNA_KIRP),
#nomatch = NA_integer_, incomparables = NULL)
#bin_SmiRNA_KIRP_underove23=bin_SmiRNA_KIRP[,temp]
#bin_SmiRNA_KIRP_underove23_sortsample=t(bin_SmiRNA_KIRP_underove23[lable_groupsample,])

#match for miRNA_KIRP_match
temp1=match(hyper_miRNA_underover_final23,colnames(miRNA_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
miRNA_KIRP_match_underove23_sortsample=t(miRNA_LGG_delmeancv[lable_groupsample,temp1])
save(miRNA_KIRP_match_underove23_sortsample,file='miRNA_KIRP_match_underove23_sortsample.Rdata')
write.csv(miRNA_KIRP_match_underove23_sortsample,"miRNA_KIRP_match_underove23_sortsample.csv",quote = F)

dim(miRNA_KIRP_match_underove23_sortsample)
#miRNA上下调
#length(hyper_miRNA_under_final[duplicated(hyper_miRNA_under_final)]) #0个重复
#去掉重复
#miRNA_over_nodup <- unique(miRNA_over$over)
hyper_miRNA_under_final=c(hyper_miRNA_underover_final[[1]],hyper_miRNA_underover_final[[3]]) 
hyper_miRNA_under_final1=Reduce(union,hyper_miRNA_under_final)
hyper_miRNA_over_final=c(hyper_miRNA_underover_final[[2]],hyper_miRNA_underover_final[[4]]) 
hyper_miRNA_over_final1=Reduce(union,hyper_miRNA_over_final)

length(unique(hyper_miRNA_under_final1))
length(unique(hyper_miRNA_over_final1))

miRNA_undergene_nodup<- unique(hyper_miRNA_under_final1)
miRNA_overgene_nodup<- unique(hyper_miRNA_over_final1)
miRNA_underovergene_nodup<- unique(c(miRNA_undergene_nodup,miRNA_overgene_nodup))
length(unique(miRNA_underovergene_nodup))
write.csv(miRNA_undergene_nodup,file='miRNA_undergene_nodup.csv')
write.csv(miRNA_overgene_nodup,file='miRNA_overgene_nodup.csv')
write.csv(miRNA_underovergene_nodup,file='miRNA_underovergene_nodup.csv')
###########################################
# heatmap for miRNA_KIRP_match_underove 2/3  #
###########################################
miRNA_differentgene_pheatmap <- miRNA_KIRP_match_underove23_sortsample
colnames(miRNA_differentgene_pheatmap) <- substr(colnames(miRNA_differentgene_pheatmap),start=1,stop=12)

setwd("F:\\ecpcs结果\\LGG下游\\差异热图")
tiff(filename = "miRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
library(pheatmap)
range(miRNA_differentgene_pheatmap)
#breaks
bk=seq(-1,1,length=6)
#with high/low labels
pheatmap(miRNA_differentgene_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=F,
         show_rownames =F, show_colnames = F,
         color = colorRampPalette(c("green",'black',"red"))(5),
         #color = colorRampPalette(c("green", "black", "firebrick3"))(7),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         border = F, # 去除cell边框
         fontsize = 8.5
)

dev.off()




######################################
# 4.2 mRNA differential analysis  #
#####################################
###########################################################
# step 1 kruskal.test for each mRNA   #
############################################################
median(SmRNA_KIRP)
quantile(SmRNA_KIRP) #可求最大最小、中位数及上下四分位数
quantile(SmRNA_KIRP, probs = c(0.25,0.75)) #求上下四分位数
setwd("F:\\ecpcs结果\\LGG下游\\differential gene/mRNA")

#choose all the feature-level association weights for the 14 subtypes
p_diffmRNA_all=NULL
for(i in 1:dim(mRNA_LGG_delmeancv)[2]){
  fit_kwH_mRNA_all=NULL
  fit_kwH_mRNA_all=kruskal.test(mRNA_LGG_delmeancv[,i]~group_sparse_kmeans)
  p_diffmRNA_all=c(p_diffmRNA_all,fit_kwH_mRNA_all$p.value)
}
length(which(p_diffmRNA_all<0.05))
save(p_diffmRNA_all,file='p_diffmRNA_all.Rdata')

#######FDR adjust########
#using the p.adjust function
q_adjust_mRNA_all <- p.adjust (p_diffmRNA_all, method = "fdr")
length(which(q_adjust_mRNA_all<0.05))
differ_mRNA_fdr=colnames(mRNA_LGG_delmeancv)[which(q_adjust_mRNA_all<0.05)]
differ_mRNA_fdr=as.data.frame(differ_mRNA_fdr)
write.csv(differ_mRNA_fdr,file='differ_mRNA_fdr.csv')
save(differ_mRNA_fdr,file='differ_mRNA_fdr.Rdata')




#############################################
#step 2 hypergeometric test for all mRNA   #
#############################################
#based on the scaled data,and choose the significant mRNA
bin_SmRNA_KIRP=SmRNA_KIRP[,which(q_adjust_mRNA_all<0.05)]
for(i in 1:dim(bin_SmRNA_KIRP)[2]){
  for(j in 1:dim(bin_SmRNA_KIRP)[1]){
    if(bin_SmRNA_KIRP[j,i]<=-0.5){
      bin_SmRNA_KIRP[j,i]=-1
    }else if(bin_SmRNA_KIRP[j,i]>=0.5){
      bin_SmRNA_KIRP[j,i]=1
    }else{
      bin_SmRNA_KIRP[j,i]=0
    }
  }
}

save(bin_SmRNA_KIRP,file='bin_SmRNA_KIRP.Rdata')

hyper_mRNA_underover_final=list()
for (m in 1:2) {
  # under expressed genes z<=-1 for cluster
  location_cluster=get(paste('location_cluster',m,sep=''))
  #location_cluster=location_cluster3
  p_hyper_cluster_under=NULL
  for(i in 1:dim(bin_SmRNA_KIRP)[2]){
    
    k_cluster=length(which(bin_SmRNA_KIRP[location_cluster,i]==-1))
    m_cluster=length(which(bin_SmRNA_KIRP[,i]==-1))
    N_cluster=dim(bin_SmRNA_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_under=c(p_hyper_cluster_under,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_under_fdr =p.adjust(p_hyper_cluster_under, method = "fdr")
  p_hyper_cluster_under_fdr_sort=sort(p_hyper_cluster_under_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_under_fdr_sortp=p_hyper_cluster_under_fdr_sort$x
  hyper_fdr_cluster_mRNA_under=colnames(bin_SmRNA_KIRP)[p_hyper_cluster_under_fdr_sort$ix][which(p_hyper_cluster_under_fdr_sortp<0.01)]
  hyper_fdr_cluster_mRNA_under=as.data.frame(hyper_fdr_cluster_mRNA_under)
  
  #for over expressed genes z>=1
  ##############################################
  # over expressed genes z>=1 for cluster  #
  ##############################################
  p_hyper_cluster_over=NULL
  for(i in 1:dim(bin_SmRNA_KIRP)[2]){
    k_cluster=length(which(bin_SmRNA_KIRP[location_cluster,i]==1))
    m_cluster=length(which(bin_SmRNA_KIRP[,i]==1))
    N_cluster=dim(bin_SmRNA_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_over=c(p_hyper_cluster_over,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_over_fdr= p.adjust(p_hyper_cluster_over, method = "fdr")
  p_hyper_cluster_over_fdr_sort=sort(p_hyper_cluster_over_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_over_fdr_sortp=p_hyper_cluster_over_fdr_sort$x
  hyper_fdr_cluster_mRNA_over=colnames(bin_SmRNA_KIRP)[p_hyper_cluster_over_fdr_sort$ix][which(p_hyper_cluster_over_fdr_sortp <0.01)]
  hyper_fdr_cluster_mRNA_over=as.data.frame(hyper_fdr_cluster_mRNA_over)
  
  
  #combine under and over
  hyper_mRNA_underover=union(hyper_fdr_cluster_mRNA_under,
                             hyper_fdr_cluster_mRNA_over)
  hyper_mRNA_underover_final=c(hyper_mRNA_underover_final,(hyper_mRNA_underover))
  
}

save(hyper_mRNA_underover_final,file='hyper_mRNA_underover_final.Rdata')
#############combine all the under and over expression genes #############
hyper_mRNA_underover_final23=Reduce(union,hyper_mRNA_underover_final)
save(hyper_mRNA_underover_final23,file='hyper_mRNA_underover_final23.Rdata')
write.table(hyper_mRNA_underover_final23,file='hyper_mRNA_underover_final23.txt',row.names = F)
#temp=match(hyper_mRNA_underover_final23,colnames(bin_SmRNA_KIRP),
#nomatch = NA_integer_, incomparables = NULL)
#bin_SmRNA_KIRP_underove23=bin_SmRNA_KIRP[,temp]
#bin_SmRNA_KIRP_underove23_sortsample=t(bin_SmRNA_KIRP_underove23[lable_groupsample,])

#match with the mRNA_KIRP_match
temp1=match(hyper_mRNA_underover_final23,colnames(mRNA_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
mRNA_KIRP_match_underove23_sortsample=t(mRNA_LGG_delmeancv[lable_groupsample,temp1])
#all(row.names(mRNA_KIRP_match_underove23_sortsample)==row.names(bin_SmRNA_KIRP_underove23_sortsample))

#save(bin_SmRNA_KIRP_underove23,file='bin_SmRNA_KIRP_underove23.Rdata')
#save(bin_SmRNA_KIRP_underove23_sortsample,file='bin_SmRNA_KIRP_underove23_sortsample.Rdata')
save(mRNA_KIRP_match_underove23_sortsample,file='mRNA_KIRP_match_underove23_sortsample.Rdata')
write.csv(mRNA_KIRP_match_underove23_sortsample,"mRNA_KIRP_match_underove23_sortsample.csv",quote = F)

dim(mRNA_KIRP_match_underove23_sortsample)

#mRNA上下调
hyper_mRNA_under_final=c(hyper_mRNA_underover_final[[1]],hyper_mRNA_underover_final[[3]]) 
hyper_mRNA_under_final1=Reduce(union,hyper_mRNA_under_final)
hyper_mRNA_over_final=c(hyper_mRNA_underover_final[[2]],hyper_mRNA_underover_final[[4]]) 
hyper_mRNA_over_final1=Reduce(union,hyper_mRNA_over_final)

length(unique(hyper_mRNA_under_final1))
length(unique(hyper_mRNA_over_final1))

mRNA_undergene_nodup<- unique(hyper_mRNA_under_final1)
mRNA_overgene_nodup<- unique(hyper_mRNA_over_final1)
mRNA_underovergene_nodup<- unique(c(mRNA_undergene_nodup,mRNA_overgene_nodup))
length(unique(mRNA_underovergene_nodup))
write.csv(mRNA_undergene_nodup,file='mRNA_undergene_nodup.csv')
write.csv(mRNA_overgene_nodup,file='mRNA_overgene_nodup.csv')
write.csv(mRNA_underovergene_nodup,file='mRNA_underovergene_nodup.csv')
###########################################
# heatmap for mRNA_KIRP_match_underove23  #
###########################################
setwd("F:\\ecpcs结果\\LGG下游\\差异热图")
mRNA_differentgene_pheatmap <- mRNA_KIRP_match_underove23_sortsample
#colnames(mRNA_differentgene_pheatmap) <- substr(colnames(mRNA_differentgene_pheatmap),start=1,stop=12)

tiff(filename = "mRNA-PRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
range(mRNA_differentgene_pheatmap)
#breaks
bk=seq(-1,1,length=6)
pheatmap(mRNA_differentgene_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=FALSE,
         show_rownames =F, show_colnames = F,
         #color = colorRampPalette(c("green",'black',"red"))(3),
         color = colorRampPalette(c("green", "black", "firebrick3"))(5),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         fontsize = 8.5,
         #cex=0.5
         #filename ='mRNA_KIRP_match_underove23_sortsample_heatmap.pdf' 
         #annotation_legend = FALSE #不显示图例
)

dev.off()




######################################
# 4.3 methy differential analysis  #
#####################################
#########################################
# step1 kruskal.test                    #
#########################################
mean(Smethy_KIRP)
quantile(Smethy_KIRP)
quantile(Smethy_KIRP, probs = c(0.25,0.75))
setwd("F:\\ecpcs结果\\LGG下游\\differential gene/methy")

p_diffmethy_all=NULL
for(i in 1:dim(methy_LGG_delmeancv)[2]){
  fit_kwH_methy_all=NULL
  fit_kwH_methy_all=kruskal.test(methy_LGG_delmeancv[,i]~group_sparse_kmeans)
  p_diffmethy_all=c(p_diffmethy_all,fit_kwH_methy_all$p.value)
}
length(which(p_diffmethy_all<0.05))
save( p_diffmethy_all,file=' p_diffmethy_all.Rdata')

#######FDR adjust########
#using the p.adjust function
q_adjust_methy_all <- p.adjust (p_diffmethy_all, method = "fdr")
length(which(q_adjust_methy_all<0.05))
differ_methy_fdr=colnames(methy_LGG_delmeancv)[which(q_adjust_methy_all<0.05)]
differ_methy_fdr=as.data.frame(differ_methy_fdr)
write.csv(differ_methy_fdr,file='differ_methy_fdr.csv')
save(differ_methy_fdr,file='differ_methy_fdr.Rdata')




###########match the related methy gene,using RefGene_Name ############
#methy_KIRP_knn_delmeancv
load("F:/TCGA_Biolink 曹老师/KIRC√/cleaned data/methy_tumor_KIRC30knn.Rdata")
methy_tumor_KIRP30knn=methy_tumor_KIRC30knn
#load("E:/English Article/Data/KIRP/methy/methy_tumor_KIRP30knn.Rdata")
location_methy_fdr=match(differ_methy_fdr[,1],methy_tumor_KIRP30knn[,1],nomatch = NA_integer_, incomparables = NULL)
methy_gene_fdr=methy_tumor_KIRP30knn[location_methy_fdr,1:4]
length(which(methy_gene_fdr$Gene_Symbol==methy_gene_fdr$RefGene_Name))
write.csv(methy_gene_fdr,file='methy_gene_fdr.csv')
save(methy_gene_fdr,file='methy_gene_fdr.Rdata')

############################################################
#step 2 hypergeometric test for all methy,cutoff=0.2,0.8   #
############################################################
#based on the scaled data,and choose the significant methylation
bin_methy_KIRP=Smethy_KIRP[,which(q_adjust_methy_all<0.05)]
for(i in 1:dim(bin_methy_KIRP)[2]){
  for(j in 1:dim(bin_methy_KIRP)[1]){
    if(bin_methy_KIRP[j,i]<=-0.5){
      bin_methy_KIRP[j,i]=-1
    }else if(bin_methy_KIRP[j,i]>=0.5){
      bin_methy_KIRP[j,i]=1
    }else{
      bin_methy_KIRP[j,i]=0
    }
  }
}
save(bin_methy_KIRP,file='bin_methy_KIRP.Rdata')

hyper_methy_underover_final=list()
for (m in 1:2) {
  # under expressed genes z<=-1 for cluster
  location_cluster=get(paste('location_cluster',m,sep=''))
  #location_cluster=location_cluster3
  p_hyper_cluster_under=NULL
  for(i in 1:dim(bin_methy_KIRP)[2]){
    
    k_cluster=length(which(bin_methy_KIRP[location_cluster,i]==-1))
    m_cluster=length(which(bin_methy_KIRP[,i]==-1))
    N_cluster=dim(bin_methy_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_under=c(p_hyper_cluster_under,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_under_fdr =p.adjust(p_hyper_cluster_under, method = "fdr")
  p_hyper_cluster_under_fdr_sort=sort(p_hyper_cluster_under_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_under_fdr_sortp=p_hyper_cluster_under_fdr_sort$x
  hyper_fdr_cluster_methy_under=colnames(bin_methy_KIRP)[p_hyper_cluster_under_fdr_sort$ix][which(p_hyper_cluster_under_fdr_sortp<0.01)]
  hyper_fdr_cluster_methy_under=as.data.frame(hyper_fdr_cluster_methy_under)
  
  #for over expressed genes z>=1
  ##############################################
  # over expressed genes z>=1 for cluster  #
  ##############################################
  p_hyper_cluster_over=NULL
  for(i in 1:dim(bin_methy_KIRP)[2]){
    k_cluster=length(which(bin_methy_KIRP[location_cluster,i]==1))
    m_cluster=length(which(bin_methy_KIRP[,i]==1))
    N_cluster=dim(bin_methy_KIRP)[1]
    n_cluster=length(location_cluster)
    p_hyper_cluster=1-phyper(k_cluster-1,m_cluster,N_cluster-m_cluster,n_cluster)
    p_hyper_cluster_over=c(p_hyper_cluster_over,p_hyper_cluster)
  }
  
  #######FDR adjust########
  #using the p.adjust function
  p_hyper_cluster_over_fdr= p.adjust(p_hyper_cluster_over, method = "fdr")
  p_hyper_cluster_over_fdr_sort=sort(p_hyper_cluster_over_fdr,decreasing = FALSE,index.return= TRUE )
  p_hyper_cluster_over_fdr_sortp= p_hyper_cluster_over_fdr_sort$x
  hyper_fdr_cluster_methy_over=colnames(bin_methy_KIRP)[p_hyper_cluster_over_fdr_sort$ix][which(p_hyper_cluster_over_fdr_sortp <0.01)]
  hyper_fdr_cluster_methy_over=as.data.frame(hyper_fdr_cluster_methy_over)
  
  
  
  #combine under and over
  hyper_methy_underover=union(hyper_fdr_cluster_methy_under,
                              hyper_fdr_cluster_methy_over)
  hyper_methy_underover_final=c(hyper_methy_underover_final,(hyper_methy_underover))
}

save(hyper_methy_underover_final,file='hyper_methy_underover_final.Rdata')

#############combine all the under and over expression genes #############
hyper_methy_underover_final23=Reduce(union,hyper_methy_underover_final)
length(hyper_methy_underover_final23)

save(hyper_methy_underover_final23,file='hyper_methy_underover_final23.Rdata')


#match with the methy_KIRP_match
temp1=match(hyper_methy_underover_final23,colnames(methy_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
methy_KIRP_match_underove23_sortsample=t(methy_LGG_delmeancv[lable_groupsample,temp1])


save(methy_KIRP_match_underove23_sortsample,file='methy_KIRP_match_underove23_sortsample.Rdata')
write.csv(methy_KIRP_match_underove23_sortsample,"methy_KIRP_match_underove23_sortsample.csv",quote = F)

dim(methy_KIRP_match_underove23_sortsample)


#甲基化上下调
hyper_methy_under_final=c(hyper_methy_underover_final[[1]],hyper_methy_underover_final[[3]]) 
hyper_methy_under_final1=Reduce(union,hyper_methy_under_final)#under去重
hyper_methy_over_final=c(hyper_methy_underover_final[[2]],hyper_methy_underover_final[[4]]) 
hyper_methy_over_final1=Reduce(union,hyper_methy_over_final)#over去重

length(unique(hyper_methy_under_final1))
length(unique(hyper_methy_over_final1))

methy_undergene_nodup<- hyper_methy_under_final1
methy_overgene_nodup<- hyper_methy_over_final1
methy_underovergene_nodup<- unique(c(methy_undergene_nodup,methy_overgene_nodup))
length(unique(methy_underovergene_nodup))
write.csv(methy_undergene_nodup,file='methy_undergene_nodup.csv')
write.csv(methy_overgene_nodup,file='methy_overgene_nodup.csv')
write.csv(methy_underovergene_nodup,file='methy_underovergene_nodup.csv')


########################################################
# heatmap for methy_KIRP_match_underove2/3_sortsample  #
#######################################################
methy_differentgene_pheatmap <- methy_KIRP_match_underove23_sortsample
methy_differentgene_pheatmap=as.matrix(methy_differentgene_pheatmap)
setwd("F:\\ecpcs结果\\LGG下游\\差异热图")
tiff(filename = "methy-PRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
range(methy_differentgene_pheatmap)
#breaks

bk=seq(-1,1,length=6)
pheatmap(methy_differentgene_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=FALSE,
         show_rownames =F, show_colnames = F,
         #color = colorRampPalette(c("green",'black',"red"))(3),
         color = colorRampPalette(c("green", "black", "firebrick3"))(5),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         fontsize = 8.5,
         border = F  # 去除cell边框
         #cex=0.5
         #filename ='mRNA_KIRP_match_underove23_sortsample_heatmap.pdf' 
)

dev.off()

###########################################
# 对差异基因进行排序，选择p最小的前15个基因  #
###########################################

###对miRNA基因进行排序，按照p值由小到大筛选前15个
setwd("F:\\ecpcs结果\\LGG下游\\top15/miRNA")
sig_gene_miRNA <- cbind(colnames(miRNA_LGG_delmeancv),q_adjust_miRNA_all)
sig_gene_miRNA <- as.data.frame(sig_gene_miRNA)
#sig_gene_sort <- as.data.frame(sig_gene[order(sig_gene[,2]),])


sig_miRNA_location <- match(rownames(miRNA_KIRP_match_underove23_sortsample),sig_gene_miRNA[,1])
sig_gene_miRNA_overlap <- as.data.frame(sig_gene_miRNA[sig_miRNA_location,])
sig_gene_miRNA_sort <- as.data.frame(sig_gene_miRNA_overlap[order(as.numeric(sig_gene_miRNA_overlap[,2]),decreasing=F),])
sig_gene_miRNA_sort_15 <- sig_gene_miRNA_sort[1:15,]

write.csv(sig_gene_miRNA_sort,file='sig_gene_miRNA_sort.csv')
write.csv(sig_gene_miRNA_sort_15,file='sig_gene_miRNA_sort_15.csv')
save(sig_gene_miRNA_sort_15,file='sig_gene_miRNA_sort_15.Rdata')

#前15热图

heat_miRNA_15_location<- match(sig_gene_miRNA_sort_15[,1],rownames(miRNA_differentgene_pheatmap))
heat_miRNA_15<-miRNA_differentgene_pheatmap[heat_miRNA_15_location,]

miRNA_differentgene15_pheatmap<-heat_miRNA_15
write.csv(miRNA_differentgene15_pheatmap,file='miRNA_differentgene15_pheatmap.csv')
miRNA_differentgene15_pheatmap=read.csv("F:\\ecpcs结果\\LGG下游\\top15\\miRNA\\miRNA_differentgene15_pheatmap.csv",header=T,stringsAsFactors = F,row.names = 1)
library(pheatmap)
tiff(filename = "miRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
range(miRNA_differentgene15_pheatmap)
#breaks
bk=seq(-1,1,length=6)
#with high/low labels
pheatmap(miRNA_differentgene15_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=F,
         show_rownames =T, show_colnames = F,
         color = colorRampPalette(c("green",'black',"red"))(6),
         #color = colorRampPalette(c("green", "black", "firebrick3"))(7),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         fontsize = 6,
         border = F,  # 去除cell边框
         annotation_legend = TRUE #不显示图例
)

dev.off()

#前15基因 kmeans  miR

#1
hsa.miR.590.3p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.590.3p"])
colnames(hsa.miR.590.3p) <- "hsa.miR.590.3p"
group_hsa.miR.590.3p <- kmeans(hsa.miR.590.3p,2,nstart = 30)
group_miR1 <- group_hsa.miR.590.3p$cluster

hsa.miR.142.3p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.142.3p"])
colnames(hsa.miR.142.3p) <- "hsa.miR.142.3p"
group_hsa.miR.142.3p <- kmeans(hsa.miR.142.3p,2,nstart = 30)
group_miR2 <- group_hsa.miR.142.3p$cluster

hsa.miR.374b <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.374b"])
colnames(hsa.miR.374b) <- "hsa.miR.374b"
group_hsa.miR.374b <- kmeans(hsa.miR.374b,2,nstart = 30)
group_miR3 <- group_hsa.miR.374b$cluster

hsa.miR.374a <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.374a"])
colnames(hsa.miR.374a) <- "hsa.miR.374a"
group_hsa.miR.374a <- kmeans(hsa.miR.374a,2,nstart = 30)
group_miR4 <- group_hsa.miR.374a$cluster

hsa.miR.675 <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.675"])
colnames(hsa.miR.675) <- "hsa.miR.675"
group_hsa.miR.675 <- kmeans(hsa.miR.675,2,nstart = 30)
group_miR5 <- group_hsa.miR.675$cluster

hsa.miR.139.3p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.139.3p"])
colnames(hsa.miR.139.3p) <- "hsa.miR.139.3p"
group_hsa.miR.139.3p <- kmeans(hsa.miR.139.3p,2,nstart = 30)
group_miR6 <- group_hsa.miR.139.3p$cluster

hsa.miR.30c.1. <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.30c.1."])
colnames(hsa.miR.30c.1.) <- "hsa.miR.30c.1."
group_hsa.miR.30c.1. <- kmeans(hsa.miR.30c.1.,2,nstart = 30)
group_miR7 <- group_hsa.miR.30c.1.$cluster

hsa.miR.139.5p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.139.5p"])
colnames(hsa.miR.139.5p) <- "hsa.miR.139.5p"
group_hsa.miR.139.5p <- kmeans(hsa.miR.139.5p,2,nstart = 30)
group_miR8 <- group_hsa.miR.139.5p$cluster

hsa.miR.328 <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.328"])
colnames(hsa.miR.328) <- "hsa.miR.328"
group_hsa.miR.328 <- kmeans(hsa.miR.328,2,nstart = 30)
group_miR9 <- group_hsa.miR.328$cluster

hsa.miR.142.5p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.142.5p"])
colnames(hsa.miR.142.5p) <- "hsa.miR.142.5p"
group_hsa.miR.142.5p <- kmeans(hsa.miR.142.5p,2,nstart = 30)
group_miR10 <- group_hsa.miR.142.5p$cluster

hsa.let.7a. <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.let.7a."])
colnames(hsa.let.7a.) <- "hsa.let.7a."
group_hsa.let.7a. <- kmeans(hsa.let.7a.,2,nstart = 30)
group_miR11 <- group_hsa.let.7a.$cluster

hsa.miR.449a <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.449a"])
colnames(hsa.miR.449a) <- "hsa.miR.449a"
group_hsa.miR.449a <- kmeans(hsa.miR.449a,2,nstart = 30)
group_miR12 <- group_hsa.miR.449a$cluster

hsa.miR.129.5p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.129.5p"])
colnames(hsa.miR.129.5p) <- "hsa.miR.129.5p"
group_hsa.miR.129.5p <- kmeans(hsa.miR.129.5p,2,nstart = 30)
group_miR13 <- group_hsa.miR.129.5p$cluster

hsa.miR.542.3p <- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.542.3p"])
colnames(hsa.miR.542.3p) <- "hsa.miR.542.3p"
group_hsa.miR.542.3p <- kmeans(hsa.miR.542.3p,2,nstart = 30)
group_miR14 <- group_hsa.miR.542.3p$cluster

hsa.miR.493.<- as.data.frame(miRNA_LGG_delmeancv[,"hsa.miR.493."])
colnames(hsa.miR.493.) <- "hsa.miR.493."
group_hsa.miR.493.<- kmeans(hsa.miR.493.,2,nstart = 30)
group_miR15 <- group_hsa.miR.493.$cluster

KIRC_miR_lable<-list()
KIRC_miR_lable[[1]]=group_miR1
KIRC_miR_lable[[2]]=group_miR2
KIRC_miR_lable[[3]]=group_miR3
KIRC_miR_lable[[4]]=group_miR4
KIRC_miR_lable[[5]]=group_miR5
KIRC_miR_lable[[6]]=group_miR6
KIRC_miR_lable[[7]]=group_miR7
KIRC_miR_lable[[8]]=group_miR8
KIRC_miR_lable[[9]]=group_miR9
KIRC_miR_lable[[10]]=group_miR10
KIRC_miR_lable[[11]]=group_miR11
KIRC_miR_lable[[12]]=group_miR12
KIRC_miR_lable[[13]]=group_miR13
KIRC_miR_lable[[14]]=group_miR14
KIRC_miR_lable[[15]]=group_miR15
p_miR =matrix(NA,15)
temp=0
for(i in 1:15){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_miR_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_miR_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_miR[[temp]]= p_value
} 
write.csv(p_miR,file='p_miR.csv')

#选123做生存曲线
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_miR1,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_miR[[1]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/miRNA")
tiff(filename = "hsa.miR.590.3p.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_miR2,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_miR[[2]],3) #3为设置的保留有效数字的位数
tiff(filename = "hsa.miR.142.3p.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_miR3,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_miR[[3]],3) #3为设置的保留有效数字的位数
tiff(filename = "hsa.miR.374b.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()


###对mRNA基因进行排序，按照p值由小到大筛选前15个
setwd("F:\\ecpcs结果\\LGG下游\\top15/mRNA")
sig_gene <- cbind(colnames(mRNA_LGG_delmeancv),q_adjust_mRNA_all)
sig_gene <- as.data.frame(sig_gene)
#sig_gene_sort <- as.data.frame(sig_gene[order(sig_gene[,2]),])


sig_overlap_location <- match(rownames(mRNA_KIRP_match_underove23_sortsample),sig_gene[,1])
sig_gene_overlap <- as.data.frame(sig_gene[sig_overlap_location,])
sig_gene_sort <- as.data.frame(sig_gene_overlap[order(sig_gene_overlap[,2]),])
sig_gene_mRNA_sort_15 <- sig_gene_sort[1:15,]

write.csv(sig_gene_mRNA_sort_15,file='sig_gene_mRNA_sort_15.csv')

#前15热图
heat_mRNA_15_location<- match(sig_gene_mRNA_sort_15[,1],rownames(mRNA_differentgene_pheatmap))
heat_mRNA_15<-mRNA_differentgene_pheatmap[heat_mRNA_15_location,]

mRNA_differentgene15_pheatmap<-heat_mRNA_15
setwd("F:\\ecpcs结果\\LGG下游\\top15/mRNA")
tiff(filename = "mRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
range(mRNA_differentgene15_pheatmap)
#breaks
bk=seq(-1,1,length=6)
#with high/low labels
pheatmap(mRNA_differentgene15_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=F,
         show_rownames =T, show_colnames = F,
         color = colorRampPalette(c("green",'black',"red"))(6),
         #color = colorRampPalette(c("green", "black", "firebrick3"))(7),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         fontsize = 6,
         border = F,  # 去除cell边框
         annotation_legend = TRUE #不显示图例
)

dev.off()

#前15基因 kmeans  mR

#1
PDZRN4<- as.data.frame(mRNA_LGG_delmeancv[,"PDZRN4"])
colnames(PDZRN4) <- "PDZRN4"
group_PDZRN4<- kmeans(PDZRN4,2,nstart = 30)
group_mR1 <- group_PDZRN4$cluster

CTAGE9<- as.data.frame(mRNA_LGG_delmeancv[,"CTAGE9"])
colnames(CTAGE9) <- "CTAGE9"
group_CTAGE9<- kmeans(CTAGE9,2,nstart = 30)
group_mR2 <- group_CTAGE9$cluster

LOC153328<- as.data.frame(mRNA_LGG_delmeancv[,"LOC153328"])
colnames(LOC153328) <- "LOC153328"
group_LOC153328<- kmeans(LOC153328,2,nstart = 30)
group_mR3 <- group_LOC153328$cluster

MSTP9<- as.data.frame(mRNA_LGG_delmeancv[,"MSTP9"])
colnames(MSTP9) <- "MSTP9"
group_MSTP9<- kmeans(MSTP9,2,nstart = 30)
group_mR4 <- group_MSTP9$cluster

C7orf41<- as.data.frame(mRNA_LGG_delmeancv[,"C7orf41"])
colnames(C7orf41) <- "C7orf41"
group_C7orf41<- kmeans(C7orf41,2,nstart = 30)
group_mR5 <- group_C7orf41$cluster

BCCIP<- as.data.frame(mRNA_LGG_delmeancv[,"BCCIP"])
colnames(BCCIP) <- "BCCIP"
group_BCCIP<- kmeans(BCCIP,2,nstart = 30)
group_mR6 <- group_BCCIP$cluster

NCOA1<- as.data.frame(mRNA_LGG_delmeancv[,"NCOA1"])
colnames(NCOA1) <- "NCOA1"
group_NCOA1<- kmeans(NCOA1,2,nstart = 30)
group_mR7 <- group_NCOA1$cluster

PRKAB2<- as.data.frame(mRNA_LGG_delmeancv[,"PRKAB2"])
colnames(PRKAB2) <- "PRKAB2"
group_PRKAB2<- kmeans(PRKAB2,2,nstart = 30)
group_mR8 <- group_PRKAB2$cluster

C3orf35<- as.data.frame(mRNA_LGG_delmeancv[,"C3orf35"])
colnames(C3orf35) <- "C3orf35"
group_C3orf35<- kmeans(C3orf35,2,nstart = 30)
group_mR9 <- group_C3orf35$cluster

TMEM87A<- as.data.frame(mRNA_LGG_delmeancv[,"TMEM87A"])
colnames(TMEM87A) <- "TMEM87A"
group_TMEM87A<- kmeans(TMEM87A,2,nstart = 30)
group_mR10 <- group_TMEM87A$cluster

C14orf102<- as.data.frame(mRNA_LGG_delmeancv[,"C14orf102"])
colnames(C14orf102) <- "C14orf102"
group_C14orf102<- kmeans(C14orf102,2,nstart = 30)
group_mR11 <- group_C14orf102$cluster

CCDC49<- as.data.frame(mRNA_LGG_delmeancv[,"CCDC49"])
colnames(CCDC49) <- "CCDC49"
group_CCDC49<- kmeans(CCDC49,2,nstart = 30)
group_mR12 <- group_CCDC49$cluster

RABGAP1L<- as.data.frame(mRNA_LGG_delmeancv[,"RABGAP1L"])
colnames(RABGAP1L) <- "RABGAP1L"
group_RABGAP1L<- kmeans(RABGAP1L,2,nstart = 30)
group_mR13 <- group_RABGAP1L$cluster

C2orf84<- as.data.frame(mRNA_LGG_delmeancv[,"C2orf84"])
colnames(C2orf84) <- "C2orf84"
group_C2orf84<- kmeans(C2orf84,2,nstart = 30)
group_mR14 <- group_C2orf84$cluster

ANKIB1<- as.data.frame(mRNA_LGG_delmeancv[,"ANKIB1"])
colnames(ANKIB1) <- "ANKIB1"
group_ANKIB1<- kmeans(ANKIB1,2,nstart = 30)
group_mR15 <- group_ANKIB1$cluster

KIRC_mR_lable<-list()
KIRC_mR_lable[[1]]=group_mR1
KIRC_mR_lable[[2]]=group_mR2
KIRC_mR_lable[[3]]=group_mR3
KIRC_mR_lable[[4]]=group_mR4
KIRC_mR_lable[[5]]=group_mR5
KIRC_mR_lable[[6]]=group_mR6
KIRC_mR_lable[[7]]=group_mR7
KIRC_mR_lable[[8]]=group_mR8
KIRC_mR_lable[[9]]=group_mR9
KIRC_mR_lable[[10]]=group_mR10
KIRC_mR_lable[[11]]=group_mR11
KIRC_mR_lable[[12]]=group_mR12
KIRC_mR_lable[[13]]=group_mR13
KIRC_mR_lable[[14]]=group_mR14
KIRC_mR_lable[[15]]=group_mR15
p_mR =matrix(NA,15)
temp=0
for(i in 1:15){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_mR_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_mR_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_mR[[temp]]= p_value
} 
write.csv(p_mR,file='p_mR.csv')

#选择6 11 12
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_mR6,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_mR[[6]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/mRNA")
tiff(filename = "BCCIP.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_mR11,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_mR[[11]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/mRNA")
tiff(filename = "C14orf102.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_mR12,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_mR[[12]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/mRNA")
tiff(filename = "CCDC49.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()


###对差异甲基化基因进行排序，选择p最小的前15个基因
setwd("F:\\ecpcs结果\\LGG下游\\top15/methy")
sig_methy <- cbind(colnames(methy_LGG_delmeancv),q_adjust_methy_all)
sig_methy <- as.data.frame(sig_methy)
#sig_gene_sort <- as.data.frame(sig_gene[order(sig_gene[,2]),])


sig_overlap_location2 <- match(rownames(methy_KIRP_match_underove23_sortsample),sig_methy[,1])
sig_methy_overlap <- as.data.frame(sig_methy[sig_overlap_location2,])
sig_methy_sort <- as.data.frame(sig_methy_overlap[order(sig_methy_overlap[,2]),])
sig_gene_methy_sort_15 <- sig_methy_sort[1:15,]
write.csv(sig_gene_methy_sort_15,file='sig_gene_methy_sort_15.csv')
#save(sig_methy_sort,file='sig_methy_sort.Rdata')

#前15热图
#热图

heat_methy_15_location<- match(sig_gene_methy_sort_15[,1],rownames(methy_differentgene_pheatmap))
heat_methy_15=methy_differentgene_pheatmap[heat_methy_15_location,]

methy_differentgene15_pheatmap<-heat_methy_15
setwd("F:\\ecpcs结果\\LGG下游\\top15/methy")
tiff(filename = "methy-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
range(methy_differentgene15_pheatmap)
#breaks
bk=seq(-1,1,length=6)
#with high/low labels
pheatmap(methy_differentgene15_pheatmap,scale='row',
         cluster_rows=T,cluster_cols=FALSE,
         show_rownames =T, show_colnames = F,
         color = colorRampPalette(c("green",'black',"red"))(6),
         #color = colorRampPalette(c("green", "black", "firebrick3"))(7),
         breaks=bk,
         legend_breaks=c(-3,3),
         legend_labels = c("Low","High"),
         annotation_col=group_sample,
         annotation_colors=colors,
         fontsize = 6,
         legend = TRUE,
         border = F,  # 去除cell边框
         annotation_legend = TRUE #不显示图例
)

dev.off()

#kmeans
TBC1D7<- as.data.frame(methy_LGG_delmeancv[,"TBC1D7"])
colnames(TBC1D7) <- "TBC1D7"
group_TBC1D7<- kmeans(TBC1D7,2,nstart = 30)
group_methy1 <- group_TBC1D7$cluster

ANP32C<- as.data.frame(methy_LGG_delmeancv[,"ANP32C"])
colnames(ANP32C) <- "ANP32C"
group_ANP32C<- kmeans(ANP32C,2,nstart = 30)
group_methy2 <- group_ANP32C$cluster

BUB1<- as.data.frame(methy_LGG_delmeancv[,"BUB1"])
colnames(BUB1) <- "BUB1"
group_BUB1<- kmeans(BUB1,2,nstart = 30)
group_methy3 <- group_BUB1$cluster

PRG4<- as.data.frame(methy_LGG_delmeancv[,"PRG4"])
colnames(PRG4) <- "PRG4"
group_PRG4<- kmeans(PRG4,2,nstart = 30)
group_methy4 <- group_PRG4$cluster

TERF2IP<- as.data.frame(methy_LGG_delmeancv[,"TERF2IP"])
colnames(TERF2IP) <- "TERF2IP"
group_TERF2IP<- kmeans(TERF2IP,2,nstart = 30)
group_methy5 <- group_TERF2IP$cluster

MTHFR<- as.data.frame(methy_LGG_delmeancv[,"MTHFR"])
colnames(MTHFR) <- "MTHFR"
group_MTHFR<- kmeans(MTHFR,2,nstart = 30)
group_methy6 <- group_MTHFR$cluster

ZBTB8OS<- as.data.frame(methy_LGG_delmeancv[,"ZBTB8OS"])
colnames(ZBTB8OS) <- "ZBTB8OS"
group_ZBTB8OS<- kmeans(ZBTB8OS,2,nstart = 30)
group_methy7 <- group_ZBTB8OS$cluster

AP4S1<- as.data.frame(methy_LGG_delmeancv[,"AP4S1"])
colnames(AP4S1) <- "AP4S1"
group_AP4S1<- kmeans(AP4S1,2,nstart = 30)
group_methy8 <- group_AP4S1$cluster

SUMO2<- as.data.frame(methy_LGG_delmeancv[,"SUMO2"])
colnames(SUMO2) <- "SUMO2"
group_SUMO2<- kmeans(SUMO2,2,nstart = 30)
group_methy9 <- group_SUMO2$cluster

NOS2A<- as.data.frame(methy_LGG_delmeancv[,"NOS2A"])
colnames(NOS2A) <- "NOS2A"
group_NOS2A<- kmeans(NOS2A,2,nstart = 30)
group_methy10 <- group_NOS2A$cluster

ISG20L1<- as.data.frame(methy_LGG_delmeancv[,"ISG20L1"])
colnames(ISG20L1) <- "ISG20L1"
group_ISG20L1<- kmeans(ISG20L1,2,nstart = 30)
group_methy11 <- group_ISG20L1$cluster

FAM101A<- as.data.frame(methy_LGG_delmeancv[,"FAM101A"])
colnames(FAM101A) <- "FAM101A"
group_FAM101A<- kmeans(FAM101A,2,nstart = 30)
group_methy12 <- group_FAM101A$cluster

POLR2I<- as.data.frame(methy_LGG_delmeancv[,"POLR2I"])
colnames(POLR2I) <- "POLR2I"
group_POLR2I<- kmeans(POLR2I,2,nstart = 30)
group_methy13 <- group_POLR2I$cluster

SLC5A5<- as.data.frame(methy_LGG_delmeancv[,"SLC5A5"])
colnames(SLC5A5) <- "SLC5A5"
group_SLC5A5<- kmeans(SLC5A5,2,nstart = 30)
group_methy14 <- group_SLC5A5$cluster

BRRN1<- as.data.frame(methy_LGG_delmeancv[,"BRRN1"])
colnames(BRRN1) <- "BRRN1"
group_BRRN1<- kmeans(BRRN1,2,nstart = 30)
group_methy15 <- group_BRRN1$cluster

KIRC_methy_lable<-list()
KIRC_methy_lable[[1]]=group_methy1
KIRC_methy_lable[[2]]=group_methy2
KIRC_methy_lable[[3]]=group_methy3
KIRC_methy_lable[[4]]=group_methy4
KIRC_methy_lable[[5]]=group_methy5
KIRC_methy_lable[[6]]=group_methy6
KIRC_methy_lable[[7]]=group_methy7
KIRC_methy_lable[[8]]=group_methy8
KIRC_methy_lable[[9]]=group_methy9
KIRC_methy_lable[[10]]=group_methy10
KIRC_methy_lable[[11]]=group_methy11
KIRC_methy_lable[[12]]=group_methy12
KIRC_methy_lable[[13]]=group_methy13
KIRC_methy_lable[[14]]=group_methy14
KIRC_methy_lable[[15]]=group_methy15
p_methy =matrix(NA,15)
temp=0
for(i in 1:15){
  KIRC_surv<-cbind(KIRC_clinical,KIRC_methy_lable[[i]])
  #手动计算P值
  temp=temp+1
  sdf=NULL
  p_value=NULL
  KIRC_surv=as.data.frame(KIRC_surv)
  sdf=survdiff(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_methy_lable[[i]], data = KIRC_surv)
  p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  p_methy[[temp]]= p_value
} 
write.csv(p_methy,file='p_methy.csv')

#选择5 1 15
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_methy5,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_methy[[5]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/methy")
tiff(filename = "TERF2IP.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_methy1,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_methy[[1]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/methy")
tiff(filename = "TBC1D7.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_methy15,
                            data= KIRC_clinical)
p_RALGPS2 <- signif(p_methy[[15]],3) #3为设置的保留有效数字的位数

mysize = 4 #p值大小设置
setwd("F:\\ecpcs结果\\LGG下游\\top15/methy")
tiff(filename = "BRRN1.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(12, "black"),
           font.y = c(12, "black"),
           font.legend=c(10,'plain','black'), #设置标签字体
           xlab='Survival time(day)',
           ylab='Survival probability',
           font.tickslab = c(12, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 4
)
dev.off()

##############################################
#
# part 8 miRNA target genes                  #
#
###############################################
setwd("F:\\ecpcs结果\\LGG下游\\go")
load("F:/课题代码/human_miRNA_target.Rdata")
miRNA_fortarget=c('hsa-miR-590-3p',
                  'hsa-miR-142-3p',
                  'hsa-miR-374b',
                  'hsa-miR-374a',
                  'hsa-miR-675',
                  'hsa-miR-139-3p',
                  'hsa-miR-30c-1',
                  'hsa-miR-139-5p',
                  'hsa-miR-328',
                  'hsa-miR-142-5p',
                  'hsa-let-7a',
                  'hsa-miR-449a',
                  'hsa-miR-129-5p',
                  'hsa-miR-542-3p',
                  'hsa-miR-493')
miRNA_fortarget=as.data.frame(miRNA_fortarget)
colnames(miRNA_fortarget)='miRNA'
miRNA_BRCA_target=merge(miRNA_fortarget, human_miRNA_target, by = 'miRNA', all= F,sort=F)
Target_Gene=unique(miRNA_BRCA_target$`Target.Gene`)

write.csv(Target_Gene,file='Target_Gene.csv')


#go富集
setwd("F:\\ecpcs结果\\LGG下游/go")
library(ggplot2)
enrich_gene=unique(c(Target_Gene,sig_gene_methy_sort_15[,1],sig_gene_mRNA_sort_15[,1]))
write.csv(enrich_gene,"enrich_gene.csv")#前15+15+靶基因

kobas_GO=read.csv("F:\\ecpcs结果\\LGG下游\\go\\lgggo10.csv")
tiff(filename = "GO2.tiff",res = 300,width = 9.5,height = 4,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family='myFont')
ggplot(kobas_GO,aes(x=richFactor,y=Term))+
  geom_point(aes(size=Input.number,color=Corrected.P.Value))+
  scale_color_gradient(low="purple",high="red")+
  labs(x="Rich factor",y="",color="P.adjust",size="Gene number")+
  theme_bw()+
  theme(legend.text=element_text(size=12,family = "myFont"),  #修改图例字体
        legend.title=element_text(size=12,family = "myFont"), #修改图例标签字体
        axis.title.x = element_text(size=13,family = "myFont"),
        axis.title.y = element_text(size=13,family = "myFont"),
        axis.text.x = element_text(size=13.5,family = "myFont"),
        axis.text.y = element_text(size=13.5,family = "myFont"))
dev.off()





kobas_kegg=read.csv("F:\\ecpcs结果\\LGG下游\\go\\lggkegg10.csv")

tiff(filename = "kegg.tiff",res = 300,width = 9.5,height = 4,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family='myFont')
ggplot(kobas_kegg,aes(x=richFactor,y=Term))+
  geom_point(aes(size=Input.number,color=Corrected.P.Value))+
  scale_color_gradient(low="purple",high="red")+
  labs(x="Rich factor",y="",color="P.adjust",size="Gene number")+
  theme_bw()+
  theme(legend.text=element_text(size=12,family = "myFont"),  #修改图例字体
        legend.title=element_text(size=12,family = "myFont"), #修改图例标签字体
        axis.title.x = element_text(size=13,family = "myFont"),
        axis.title.y = element_text(size=13,family = "myFont"),
        axis.text.x = element_text(size=13.5,family = "myFont"),
        axis.text.y = element_text(size=13.5,family = "myFont"))
dev.off()



#######################################################
#
# pathway activities analysis                #
#
#######################################################
library(progeny)
#library(ggbiplot)
library(RColorBrewer)
library(igraph)
#########################################################
# ECPCS #
#########################################################

#########################################################
# step1 calculate pathway activities for mRNA_KIRP_match #
#########################################################
mRNA_KIRP_tmatch <- as.matrix(t(mRNA_LGG_delmeancv)) #progeny()函数要求输入数据必须为矩阵,行为基因，列为样本
score_pathways <- progeny(mRNA_KIRP_tmatch, scale=TRUE, organism="Human", top = 100, perm = 1)


#########################################################
# step2 kruskal.test for pathway activities in different clusters #
#########################################################
colnames(score_pathways)
p_pathway_all=NULL
for(i in 1:dim(score_pathways)[2]){
  fit_kwH_pathway=NULL
  fit_kwH_pathway=kruskal.test(score_pathways[,i]~as.factor(group_sparse_kmeans))
  p_pathway_all=c(p_pathway_all,fit_kwH_pathway$p.value)
}

fit_kwH_pathway
#######FDR adjust########
#using the p.adjust function
q_adjust_pathway_all <- p.adjust(p_pathway_all, method = "fdr")
length(which(q_adjust_pathway_all<0.05))

differ_pathway_fdr=as.data.frame(q_adjust_pathway_all)
row.names(differ_pathway_fdr)=colnames(score_pathways)
differ_pathway_fdr_sort=sort(differ_pathway_fdr[,1],decreasing =F,index.return=T )
differ_pathway_fdr_sortp=differ_pathway_fdr[differ_pathway_fdr_sort$ix,1,drop = FALSE]

setwd("F:\\ecpcs结果\\LGG下游/pathway activity")
save(differ_pathway_fdr,file='differ_pathway_fdr.Rdata')
write.csv(differ_pathway_fdr,file='differ_pathway_fdr.csv')


#########################################################
# step3 boxplot for average pathway activities in different clusters #
#########################################################

library(ggplot2)
library(ggpubr)

pathway_box <- cbind(score_pathways,group_sparse_kmeans)

## method 1 for boxplot
mainname=row.names(differ_pathway_fdr)[3] #可以修改g，使箱式图分别显示不同通路结果
mainname

EGFR <- as.data.frame(pathway_box[,c("EGFR","group_sparse_kmeans")])
colnames(EGFR) <- c("Pathway activity","Cluster")
VEGF <- as.data.frame(pathway_box[,c("VEGF","group_sparse_kmeans")])
colnames(VEGF) <- c("Pathway activity","Cluster")
MAPK <- as.data.frame(pathway_box[,c("MAPK","group_sparse_kmeans")])
colnames(MAPK) <- c("Pathway activity","Cluster")
Androgen <- as.data.frame(pathway_box[,c("Androgen","group_sparse_kmeans")])
colnames(Androgen) <- c("Pathway activity","Cluster")
Estrogen <- as.data.frame(pathway_box[,c("Estrogen","group_sparse_kmeans")])
colnames(Estrogen) <- c("Pathway activity","Cluster")
TNFα <- as.data.frame(pathway_box[,c("TNFa","group_sparse_kmeans")])
colnames(TNFα) <- c("Pathway activity","Cluster")
p53 <- as.data.frame(pathway_box[,c("p53","group_sparse_kmeans")])
colnames(p53) <- c("Pathway activity","Cluster")


setwd("F:\\ecpcs结果\\LGG下游/pathway activity")

mysize = 4 #p值大小设置
tiff(filename = "EGFR.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(EGFR,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "EGFR",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.9,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "VEGF.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggboxplot(VEGF,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "VEGF",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=3.8,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "MAPK.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggboxplot(MAPK,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "MAPK",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=3,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "Androgen.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggboxplot(Androgen,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "Androgen",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=1.8,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()


tiff(filename = "Estrogen.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggboxplot(Estrogen,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "Estrogen",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.0,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "TNFα.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggboxplot(TNFα,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "TNFα",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.6,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "p53.tiff",res = 300,width = 3,height = 3,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggboxplot(p53,x='Cluster',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "p53",font.title=12,
          font.x = c(12, "black"),font.y = c(12, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.6,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 10), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 10)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()
