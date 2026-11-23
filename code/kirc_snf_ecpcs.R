library(MASS)
library(SNFtool)
library("Matrix")
library("parallel")
library(survival)
library(survminer)
library("R.matlab")

#设置为当前代码路径
setwd("F:\\课题代码")

load("F:/data/老师的数据和分析/real data analysis 1st/KIRC/final cleaned data keep race3/clinic_KIRC_keeprace3.Rdata")

KIRC_clinical=clinic_KIRC_keeprace3

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

library(MOVICS)
moic1.res.list=NULL
samID=row.names(KIRC_clinical)
moic1.res.list[["SNF"]][["fit"]]=W_snf
clust.res.snf=data.frame(samID,clust=group3_snf)
moic1.res.list[["SNF"]][["clust.res"]]=clust.res.snf
moic1.res.list[["SNF"]][["mo.method"]]="SNF"

moic1.res.list[["NE"]][["fit"]]=W_NE
clust.res.NE=data.frame(samID,clust=group4_NE)
moic1.res.list[["NE"]][["clust.res"]]=clust.res.NE
moic1.res.list[["NE"]][["mo.method"]]="NE"

moic1.res.list[["skf"]][["fit"]]=W_skf
clust.res.skf=data.frame(samID,clust=group2_skf)
moic1.res.list[["skf"]][["clust.res"]]=clust.res.skf
moic1.res.list[["skf"]][["mo.method"]]="skf"

moic1.res.list[["hopes"]][["fit"]]=W_hopes
clust.res.hopes=data.frame(samID,clust=group3_hopes)
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
#################################################################################################
#
#part  生存曲线
#
##############################################################################################

setwd("F:\\ecpcs结果\\KIRC下游")
#未合并
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ KIRC_ECPCS_lable[[2]], 
             data = KIRC_surv)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_value <- signif(p_value,3)
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
tiff(filename = "KIRC1.tiff",res = 300,width = 5,height = 3.8,units = 'in',compression = c("lzw"))
fit<-survfit(Surv(as.numeric(overall_survival),as.numeric(vital_status))~KIRC_ECPCS_lable[[2]],data=KIRC_surv)
ggsurvplot(fit, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2","Cluster3"), # 设置图例
           palette = c("#FF5F5F", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_value, # 显示 pvalue
           pval.size = 6
)

dev.off()

#合并生存曲线
sdf=NULL
p_value=NULL
e2=KIRC_ECPCS_lable[[2]]
e2[which(e2=='3')]=1
group_sparse_kmeans=e2
write.csv(group_sparse_kmeans,"group_sparse_kmeans.csv")
group_gaodi=read.csv("group_gaodi.csv",header=T,stringsAsFactors = F)

group_gaodi <- as.matrix(group_gaodi)

#group_gaodi <- as.integer(group_gaodi)


sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ e2, 
             data = KIRC_surv)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_value <- signif(p_value,3)
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
tiff(filename = "KIRC.tiff",res = 300,width = 5,height = 3.8,units = 'in',compression = c("lzw"))
fit<-survfit(Surv(as.numeric(overall_survival),as.numeric(vital_status))~e2,data=KIRC_surv)
ggsurvplot(fit, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Low-risk group", "High-risk group"), # 设置图例
           palette = c("#FF5F5F", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_value, # 显示 pvalue
           pval.size = 6
)

dev.off()
#################################################################################################
#
#part5 寻找差异基因
#
################################################################################################
#load("F:/课题代码/MATLAB TCGA DATA/CHONGZHENG/中国LGG数据/quant10 3omics/zhong.RData")


group_sample_sort=sort(group_gaodi,decreasing = FALSE,index.return=T)
lable_groupsample=group_sample_sort$ix
group_sample=data.frame(as.character(group_sample_sort$x))# sample分组文件
rownames(group_sample)=rownames(KIRC_clinical)[lable_groupsample] #将排序后的分组文件行名也按照分组排列顺序重新排列
colnames(group_sample)='Group'
#write.csv(group_sample,"group_sample.csv")

#colors=list(Cluster=c('1'="#ff7473", '2'="#4bc2c5",'3'="#ffc33c"))
colors=list(Cluster=c('1'="#ff7473",'2'="#4bc2c5"))

###############################################################################
# location for each cluster #
###############################################################################
table(group_sparse_kmeans)

location_cluster1=which(group_sparse_kmeans==1)
location_cluster2=which(group_sparse_kmeans==2)


group_sparse_kmeans <- as.matrix(group_sparse_kmeans)
group_sparse_kmeans <- as.integer(group_sparse_kmeans)
## 定义下调上调表达量
count_under=function(x){
  length(which(x==-1))
}

count_over=function(x){
  length(which(x==1))
}

load("F:/老师的数据和分析/real data analysis 1st/KIRC/final cleaned data keep race3/methy_KIRC_match.Rdata")
load("F:/老师的数据和分析/real data analysis 1st/KIRC/final cleaned data keep race3/miRNA_KIRC_match.Rdata")
load("F:/老师的数据和分析/real data analysis 1st/KIRC/final cleaned data keep race3/mRNA_KIRC_match.Rdata")

miRNA_LGG_delmeancv=miRNA_KIRC_match
mRNA_LGG_delmeancv=mRNA_KIRC_match
methy_LGG_delmeancv=methy_KIRC_match
SmiRNA_KIRP=standardNormalization(miRNA_LGG_delmeancv)
SmRNA_KIRP=standardNormalization(mRNA_LGG_delmeancv)
Smethy_KIRP=standardNormalization(methy_LGG_delmeancv)
########################################################################################################
# 4.1 miRNA differential analysis  #
#
#########################################################################################################

#quantile(SmiRNA_KIRP, probs = c(0.25,0.75))
#quantile(SmiRNA_KIRP)
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/miRNA")

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
length(which(q_adjust_miRNA_all<0.01))
differ_miRNA_fdr=colnames(miRNA_LGG_delmeancv)[which(q_adjust_miRNA_all<0.01)]
differ_miRNA_fdr=as.data.frame(differ_miRNA_fdr)
#write.csv(differ_miRNA_fdr,file='differ_miRNA_fdr.csv')
save(differ_miRNA_fdr,file='differ_miRNA_fdr.Rdata')



########################################################################################
#step 2 hypergeometric test for the fdar significant miRNA   #
########################################################################################
length(which(q_adjust_miRNA_all>=0.01))

#based on the scaled data,and choose the significant miRNA 定义<0.5为-1，>=0.5为1
bin_SmiRNA_KIRP=SmiRNA_KIRP[,which(q_adjust_miRNA_all<0.01)]
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
temp=match(hyper_miRNA_underover_final23,colnames(bin_SmiRNA_KIRP),
           nomatch = NA_integer_, incomparables = NULL)
bin_SmiRNA_KIRP_underove23=bin_SmiRNA_KIRP[,temp]
bin_SmiRNA_KIRP_underove23_sortsample=t(bin_SmiRNA_KIRP_underove23[lable_groupsample,])

#match for miRNA_KIRP_match
temp1=match(hyper_miRNA_underover_final23,colnames(miRNA_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
miRNA_KIRP_match_underove23_sortsample=t(miRNA_LGG_delmeancv[lable_groupsample,temp1])

###对基因进行排序，按照p值由小到大筛选前15个
sig_gene_miRNA <- cbind(colnames(miRNA_LGG_delmeancv),q_adjust_miRNA_all)
sig_gene_miRNA <- as.data.frame(sig_gene_miRNA)
#sig_gene_sort <- as.data.frame(sig_gene[order(sig_gene[,2]),])


sig_miRNA_location <- match(rownames(miRNA_KIRP_match_underove23_sortsample),sig_gene_miRNA[,1])
sig_gene_miRNA_overlap <- as.data.frame(sig_gene_miRNA[sig_miRNA_location,])
sig_gene_miRNA_sort <- as.data.frame(sig_gene_miRNA_overlap[order(as.numeric(sig_gene_miRNA_overlap[,2]),decreasing=F),])
sig_gene_miRNA_sort_15 <- sig_gene_miRNA_sort[1:15,]

write.csv(sig_gene_miRNA_sort_15,file='sig_gene_miRNA_sort_15.csv')
save(sig_gene_miRNA_sort_15,file='sig_gene_miRNA_sort_15.Rdata')

save(bin_SmiRNA_KIRP_underove23,file='bin_SmiRNA_KIRP_underove23.Rdata')
save(bin_SmiRNA_KIRP_underove23_sortsample,file='bin_SmiRNA_KIRP_underove23_sortsample.Rdata')
save(miRNA_KIRP_match_underove23_sortsample,file='miRNA_KIRP_match_underove23_sortsample.Rdata')
write.csv(miRNA_KIRP_match_underove23_sortsample,"miRNA_KIRP_match_underove23_sortsample.csv",quote = F)

dim(miRNA_KIRP_match_underove23_sortsample)

###########################################
# heatmap for miRNA_KIRP_match_underove 2/3  #
###########################################
miRNA_differentgene_pheatmap <- miRNA_KIRP_match_underove23_sortsample
colnames(miRNA_differentgene_pheatmap) <- substr(colnames(miRNA_differentgene_pheatmap),start=1,stop=12)

setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/miRNA")
tiff(filename = "miRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
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
         annotation_legend = FALSE #不显示图例
         #fontsize = 10
)

dev.off()


setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/miRNA")
tiff(filename = "miRNA-RRCC1.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
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
         fontsize = 10
         
)

dev.off()
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
write.csv(miRNA_undergene_nodup,file='miRNA_undergene_nodup.csv')
write.csv(miRNA_overgene_nodup,file='miRNA_overgene_nodup.csv')
write.csv(miRNA_underovergene_nodup,file='miRNA_underovergene_nodup.csv')



######################################
# 4.2 mRNA differential analysis  #
#####################################
###########################################################
# step 1 kruskal.test for each mRNA   #
############################################################
median(SmRNA_KIRP)
quantile(SmRNA_KIRP) #可求最大最小、中位数及上下四分位数
quantile(SmRNA_KIRP, probs = c(0.25,0.75)) #求上下四分位数
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/mRNA")

#choose all the feature-level association weights for the 14 subtypes
p_diffmRNA_all=NULL
for(i in 1:dim(mRNA_LGG_delmeancv)[2]){
  fit_kwH_mRNA_all=NULL
  fit_kwH_mRNA_all=kruskal.test(mRNA_LGG_delmeancv[,i]~group_sparse_kmeans)
  p_diffmRNA_all=c(p_diffmRNA_all,fit_kwH_mRNA_all$p.value)
}
length(which(p_diffmRNA_all<0.01))
save(p_diffmRNA_all,file='p_diffmRNA_all.Rdata')

#######FDR adjust########
#using the p.adjust function
q_adjust_mRNA_all <- p.adjust (p_diffmRNA_all, method = "fdr")
length(which(q_adjust_mRNA_all<0.01))
differ_mRNA_fdr=colnames(mRNA_LGG_delmeancv)[which(q_adjust_mRNA_all<0.01)]
differ_mRNA_fdr=as.data.frame(differ_mRNA_fdr)
write.csv(differ_mRNA_fdr,file='differ_mRNA_fdr.csv')
save(differ_mRNA_fdr,file='differ_mRNA_fdr.Rdata')




#############################################
#step 2 hypergeometric test for all mRNA   #
#############################################
#based on the scaled data,and choose the significant mRNA
bin_SmRNA_KIRP=SmRNA_KIRP[,which(q_adjust_mRNA_all<0.01)]
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
temp=match(hyper_mRNA_underover_final23,colnames(bin_SmRNA_KIRP),
           nomatch = NA_integer_, incomparables = NULL)
bin_SmRNA_KIRP_underove23=bin_SmRNA_KIRP[,temp]
bin_SmRNA_KIRP_underove23_sortsample=t(bin_SmRNA_KIRP_underove23[lable_groupsample,])

#match with the mRNA_KIRP_match
temp1=match(hyper_mRNA_underover_final23,colnames(mRNA_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
mRNA_KIRP_match_underove23_sortsample=t(mRNA_LGG_delmeancv[lable_groupsample,temp1])
all(row.names(mRNA_KIRP_match_underove23_sortsample)==row.names(bin_SmRNA_KIRP_underove23_sortsample))

save(bin_SmRNA_KIRP_underove23,file='bin_SmRNA_KIRP_underove23.Rdata')
save(bin_SmRNA_KIRP_underove23_sortsample,file='bin_SmRNA_KIRP_underove23_sortsample.Rdata')
save(mRNA_KIRP_match_underove23_sortsample,file='mRNA_KIRP_match_underove23_sortsample.Rdata')
write.csv(mRNA_KIRP_match_underove23_sortsample,"mRNA_KIRP_match_underove23_sortsample.csv",quote = F)

dim(mRNA_KIRP_match_underove23_sortsample)

###对基因进行排序，按照p值由小到大筛选前15个
sig_gene <- cbind(colnames(mRNA_LGG_delmeancv),q_adjust_mRNA_all)
sig_gene <- as.data.frame(sig_gene)
#sig_gene_sort <- as.data.frame(sig_gene[order(sig_gene[,2]),])


sig_overlap_location <- match(rownames(mRNA_KIRP_match_underove23_sortsample),sig_gene[,1])
sig_gene_overlap <- as.data.frame(sig_gene[sig_overlap_location,])
sig_gene_mRNA_sort <- as.data.frame(sig_gene_overlap[order(as.numeric(sig_gene_overlap[,2]),decreasing=F),])
sig_gene_mRNA_sort_15 <- sig_gene_mRNA_sort[1:15,]
write.csv(sig_gene_mRNA_sort,file='sig_gene_mRNA_sort.csv')
write.csv(sig_gene_mRNA_sort_15,file='sig_gene_mRNA_sort_15.csv')
save(sig_gene_mRNA_sort,file='sig_gene_mRNA_sort.Rdata')
###########################################
# heatmap for mRNA_KIRP_match_underove23  #
###########################################
mRNA_differentgene_pheatmap <- mRNA_KIRP_match_underove23_sortsample
#colnames(mRNA_differentgene_pheatmap) <- substr(colnames(mRNA_differentgene_pheatmap),start=1,stop=12)
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/mRNA")
tiff(filename = "mRNA-PRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
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
         #fontsize = 8.5,
         #cex=0.5
         #filename ='mRNA_KIRP_match_underove23_sortsample_heatmap.pdf' 
         annotation_legend = FALSE #不显示图例
)

dev.off()

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
write.csv(methy_undergene_nodup,file='methy_undergene_nodup.csv')
write.csv(methy_overgene_nodup,file='methy_overgene_nodup.csv')
write.csv(methy_underovergene_nodup,file='methy_underovergene_nodup.csv')


######################################
# 4.3 methy differential analysis  #
#####################################
#########################################
# step1 kruskal.test                    #
#########################################
mean(Smethy_KIRP)
quantile(Smethy_KIRP)
quantile(Smethy_KIRP, probs = c(0.25,0.75))
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/methy")

p_diffmethy_all=NULL
for(i in 1:dim(methy_LGG_delmeancv)[2]){
  fit_kwH_methy_all=NULL
  fit_kwH_methy_all=kruskal.test(methy_LGG_delmeancv[,i]~group_sparse_kmeans)
  p_diffmethy_all=c(p_diffmethy_all,fit_kwH_methy_all$p.value)
}
length(which(p_diffmethy_all<0.01))
save( p_diffmethy_all,file=' p_diffmethy_all.Rdata')

#######FDR adjust########
#using the p.adjust function
q_adjust_methy_all <- p.adjust (p_diffmethy_all, method = "fdr")
length(which(q_adjust_methy_all<0.01))
differ_methy_fdr=colnames(methy_LGG_delmeancv)[which(q_adjust_methy_all<0.01)]
differ_methy_fdr=as.data.frame(differ_methy_fdr)
write.csv(differ_methy_fdr,file='differ_methy_fdr.csv')
save(differ_methy_fdr,file='differ_methy_fdr.Rdata')




###########match the related methy gene,using RefGene_Name ############
#methy_KIRP_knn_delmeancv
load("F:/TCGA_Biolink 曹老师/KIRC√/cleaned data/methy_tumor_KIRC30knn.Rdata")
methy_tumor_KIRP30knn=methy_tumor_KIRC30knn

location_methy_fdr=match(differ_methy_fdr[,1],methy_tumor_KIRP30knn[,1],nomatch = NA_integer_, incomparables = NULL)
methy_gene_fdr=methy_tumor_KIRP30knn[location_methy_fdr,1:4]
length(which(methy_gene_fdr$Gene_Symbol==methy_gene_fdr$RefGene_Name))
write.csv(methy_gene_fdr,file='methy_gene_fdr.csv')
save(methy_gene_fdr,file='methy_gene_fdr.Rdata')

############################################################
#step 2 hypergeometric test for all methy,cutoff=0.2,0.8   #
############################################################
#based on the scaled data,and choose the significant methylation
bin_methy_KIRP=Smethy_KIRP[,which(q_adjust_methy_all<0.01)]
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
length(hyper_methy_underover_final23)#总位点

save(hyper_methy_underover_final23,file='hyper_methy_underover_final23.Rdata')

temp=match(hyper_methy_underover_final23,colnames(bin_methy_KIRP),
           nomatch = NA_integer_, incomparables = NULL)
bin_methy_KIRP_underove23=bin_methy_KIRP[,temp]
bin_methy_KIRP_underove23_sortsample=t(bin_methy_KIRP_underove23[lable_groupsample,])

#match with the methy_KIRP_match
temp1=match(hyper_methy_underover_final23,colnames(methy_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
methy_KIRP_match_underove23_sortsample=t(methy_LGG_delmeancv[lable_groupsample,temp1])
all(row.names(bin_methy_KIRP_underove23_sortsample)==row.names(methy_KIRP_match_underove23_sortsample))
methy_gene_fdrenrich=methy_tumor_KIRP30knn[temp1,1:4]
methy_gene_fdrenrich_symbol=methy_tumor_KIRP30knn[temp1,3]
write.table(methy_gene_fdrenrich,file='methy_gene_fdrenrich.txt',row.names = F)
write.table(methy_gene_fdrenrich_symbol,file='methy_gene_fdrenrich_symbol.txt',row.names = F)
length(which(methy_gene_fdrenrich[,4]==''))
length(which(methy_gene_fdrenrich[,3]==''))
length(unique(methy_gene_fdrenrich[,3]))

save(bin_methy_KIRP_underove23,file='bin_methy_KIRP_underove23.Rdata')
save(bin_methy_KIRP_underove23_sortsample,file='bin_methy_KIRP_underove23_sortsample.Rdata')
save(methy_KIRP_match_underove23_sortsample,file='methy_KIRP_match_underove23_sortsample.Rdata')
write.csv(methy_KIRP_match_underove23_sortsample,"methy_KIRP_match_underove23_sortsample.csv",quote = F)

dim(methy_KIRP_match_underove23_sortsample)

########################################################
# heatmap for methy_KIRP_match_underove2/3_sortsample  #
#######################################################
methy_differentgene_pheatmap <- methy_KIRP_match_underove23_sortsample
colnames(methy_differentgene_pheatmap) <- substr(colnames(methy_differentgene_pheatmap),start=1,stop=12)
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/methy")
tiff(filename = "methy-PRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
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
         #fontsize = 8.5
         #cex=0.5
         #filename ='mRNA_KIRP_match_underove23_sortsample_heatmap.pdf' 
         annotation_legend = FALSE #不显示图例
)


dev.off()
#对差异甲基化基因进行排序，选择p最小的前15个基因
sig_methy <- cbind(methy_tumor_KIRC30knn[,c(1,3)],q_adjust_methy_all)
sig_methy <- as.data.frame(sig_methy)


sig_overlap_location2 <- match(rownames(methy_KIRP_match_underove23_sortsample),sig_methy[,1])
sig_methy_overlap <- as.data.frame(sig_methy[sig_overlap_location2,])
sig_methy_sort <- as.data.frame(sig_methy_overlap[order(as.numeric(sig_methy_overlap[,2]),decreasing=F),])
sig_methy_sort_15 <- sig_methy_sort[1:15,]


#将前15个位点对应基因

load("F:/TCGA_Biolink 曹老师/KIRC√/cleaned data/methy_tumor_KIRC30knn.Rdata")
methy_tumor_KIRP30knn=methy_tumor_KIRC30knn

location_sig <- match(sig_methy_sort_15[,1],methy_tumor_KIRP30knn[,1],nomatch = NA_integer_, incomparables = NULL) #返回差异基因在原始数据中的位置
methy_different_sig <- methy_tumor_KIRP30knn[location_sig,] #提取差异基因
sig_gene_methy_sort_15=methy_different_sig[,3]

write.csv(sig_methy_sort,file='sig_methy_sort.csv')
write.csv(sig_gene_methy_sort_15,file='sig_gene_methy_sort_15.csv')




#甲基化上下调
hyper_methy_under_final=c(hyper_methy_underover_final[[1]],hyper_methy_underover_final[[3]]) 
hyper_methy_under_final1=Reduce(union,hyper_methy_under_final)#位点under去重
hyper_methy_over_final=c(hyper_methy_underover_final[[2]],hyper_methy_underover_final[[4]]) 
hyper_methy_over_final1=Reduce(union,hyper_methy_over_final)#位点over去重
temp3=match(hyper_methy_under_final1,colnames(methy_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)#匹配基因
methy_undergene_fdrenrich_symbol=methy_tumor_KIRP30knn[temp3,3]
temp4=match(hyper_methy_over_final1,colnames(methy_LGG_delmeancv),
            nomatch = NA_integer_, incomparables = NULL)
methy_overgene_fdrenrich_symbol=methy_tumor_KIRP30knn[temp4,3]
length(unique(methy_undergene_fdrenrich_symbol))
length(unique(methy_overgene_fdrenrich_symbol))

methy_undergene_nodup<- unique(methy_undergene_fdrenrich_symbol)#基因under去重
methy_overgene_nodup<- unique(methy_overgene_fdrenrich_symbol)
methy_underovergene_nodup<- unique(c(methy_undergene_nodup,methy_overgene_nodup))
write.csv(methy_undergene_nodup,file='methy_undergene_nodup.csv')
write.csv(methy_overgene_nodup,file='methy_overgene_nodup.csv')
write.csv(methy_underovergene_nodup,file='methy_underovergene_nodup.csv')

##############################################
#
#  前15基因 kmeans
###############################################
#mRNA
#TMEM150C
setwd("F:\\ecpcs结果\\KIRC下游/top15\\mrna")
#热图

heat_mRNA_15_location<- match(sig_gene_mRNA_sort_15[,1],rownames(mRNA_differentgene_pheatmap))
heat_mRNA_15<-mRNA_differentgene_pheatmap[heat_mRNA_15_location,]

mRNA_differentgene15_pheatmap<-heat_mRNA_15

tiff(filename = "mRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
library(pheatmap)
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
         fontsize = 6.5,
         annotation_legend = TRUE #不显示图例
)

dev.off()

TMEM150C <- as.data.frame(mRNA_LGG_delmeancv[,"TMEM150C"])
colnames(TMEM150C) <- "TMEM150C"
group_RALGPS2 <- kmeans(TMEM150C,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "TMEM150C.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()


#	NPR3
NPR3 <- as.data.frame(mRNA_LGG_delmeancv[,"NPR3"])
colnames(NPR3) <- "NPR3"
group_RALGPS2 <- kmeans(NPR3,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "NPR3.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#	BCL2
BCL2 <- as.data.frame(mRNA_LGG_delmeancv[,"BCL2"])
colnames(BCL2) <- "BCL2"
group_RALGPS2 <- kmeans(BCL2,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "BCL2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#CAVIN2
CAVIN2 <- as.data.frame(mRNA_LGG_delmeancv[,"CAVIN2"])
colnames(CAVIN2) <- "CAVIN2"
group_RALGPS2 <- kmeans(CAVIN2,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "CAVIN2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#NTN4
NTN4<- as.data.frame(mRNA_LGG_delmeancv[,"NTN4"])
colnames(NTN4) <- "NTN4"
group_RALGPS2 <- kmeans(NTN4,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "NTN4.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#6
EDNRB<- as.data.frame(mRNA_LGG_delmeancv[,"EDNRB"])
colnames(EDNRB) <- "EDNRB"
group_RALGPS2 <- kmeans(EDNRB,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "EDNRB.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#7
PTPRH<- as.data.frame(mRNA_LGG_delmeancv[,"PTPRH"])
colnames(PTPRH) <- "PTPRH"
group_RALGPS2 <- kmeans(PTPRH,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "PTPRH.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#8
LIMK1<- as.data.frame(mRNA_LGG_delmeancv[,"LIMK1"])
colnames(LIMK1) <- "LIMK1"
group_RALGPS2 <- kmeans(LIMK1,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "LIMK1.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#9
PDGFD<- as.data.frame(mRNA_LGG_delmeancv[,"PDGFD"])
colnames(PDGFD) <- "PDGFD"
group_RALGPS2 <- kmeans(PDGFD,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "PDGFD.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#10
AQP1<- as.data.frame(mRNA_LGG_delmeancv[,"AQP1"])
colnames(AQP1) <- "AQP1"
group_RALGPS2 <- kmeans(AQP1,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "AQP1.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#11
PKHD1<- as.data.frame(mRNA_LGG_delmeancv[,"PKHD1"])
colnames(PKHD1) <- "PKHD1"
group_RALGPS2 <- kmeans(PKHD1,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "PKHD1.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#12
SLC27A2<- as.data.frame(mRNA_LGG_delmeancv[,"SLC27A2"])
colnames(SLC27A2) <- "SLC27A2"
group_RALGPS2 <- kmeans(SLC27A2,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "SLC27A2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#13
MYLIP<- as.data.frame(mRNA_LGG_delmeancv[,"MYLIP"])
colnames(MYLIP) <- "MYLIP"
group_RALGPS2 <- kmeans(MYLIP,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "MYLIP.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#14
PLPP3<- as.data.frame(mRNA_LGG_delmeancv[,"PLPP3"])
colnames(PLPP3) <- "PLPP3"
group_RALGPS2 <- kmeans(PLPP3,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "PLPP3.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#15
TRPM3<- as.data.frame(mRNA_LGG_delmeancv[,"TRPM3"])
colnames(TRPM3) <- "TRPM3"
group_RALGPS2 <- kmeans(TRPM3,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "TRPM3.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()


#miRNA
#  1 hsa-mir-139
setwd("F:\\ecpcs结果\\KIRC下游/top15\\mirna")
heat_miRNA_15_location<- match(sig_gene_miRNA_sort_15[,1],rownames(miRNA_differentgene_pheatmap))
heat_miRNA_15<-miRNA_differentgene_pheatmap[heat_miRNA_15_location,]

miRNA_differentgene15_pheatmap<-heat_miRNA_15

tiff(filename = "miRNA-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
library(pheatmap)
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
         fontsize = 6.5,
         annotation_legend = TRUE #不显示图例
)

dev.off()

hsa_mir_139<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-139"])
colnames(hsa_mir_139) <- "hsa_mir_139"
group_RALGPS2 <- kmeans(hsa_mir_139,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置



tiff(filename = "hsa_mir_139.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman

ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#2
hsa_mir_21<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-21"])
colnames(hsa_mir_21) <- "hsa_mir_21"
group_RALGPS2 <- kmeans(hsa_mir_21,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa_mir_21.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#3
hsa_mir_204<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-204"])
colnames(hsa_mir_204) <- "hsa_mir_204"
group_RALGPS2 <- kmeans(hsa_mir_204,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa_mir_204.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#4
hsa_mir_10<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-10b"])
colnames(hsa_mir_10) <- "hsa-mir-10b"
group_RALGPS2 <- kmeans(hsa_mir_10,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-10b.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#5
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-1468"])
colnames(hsa_mir_1468) <- "hsa-mir-1468"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-1468.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#6
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-130b"])
colnames(hsa_mir_1468) <- "hsa-mir-130b"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-130b.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#7
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-625"])
colnames(hsa_mir_1468) <- "hsa-mir-625"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-625.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#8
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-221"])
colnames(hsa_mir_1468) <- "hsa-mir-221"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-221.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#9
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-193b"])
colnames(hsa_mir_1468) <- "hsa-mir-193b"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-193b.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#10
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-9-1"])
colnames(hsa_mir_1468) <- "hsa-mir-9-1"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-9-1.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()

#11
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-9-2"])
colnames(hsa_mir_1468) <- "hsa-mir-9-2"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-9-2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#12
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-9-3"])
colnames(hsa_mir_1468) <- "hsa-mir-9-3"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-9-3.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#13
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-671"])
colnames(hsa_mir_1468) <- "hsa-mir-671"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-671.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#14
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-27b"])
colnames(hsa_mir_1468) <- "hsa-mir-27b"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-27b.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()
#15
hsa_mir_1468<- as.data.frame(miRNA_LGG_delmeancv[,"hsa-mir-365a"])
colnames(hsa_mir_1468) <- "hsa-mir-365a"
group_RALGPS2 <- kmeans(hsa_mir_1468,2,nstart = 30)
group_RALGPS2_surv <- group_RALGPS2$cluster
fit_RALGPS2_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~group_RALGPS2_surv,
                            data= KIRC_clinical)
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RALGPS2_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RALGPS2=p_value
p_RALGPS2 <- signif(p_value,3) #3为设置的保留有效数字的位数

mysize = 6 #p值大小设置

tiff(filename = "hsa-mir-365a.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RALGPS2_surv, linetype = "strata",
           #pval=T,
           #size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5","#ffc33c"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RALGPS2, # 显示 pvalue
           pval.size = 6
)
dev.off()


###### for methy
setwd("F:\\ecpcs结果\\KIRC下游/top15\\methy")

#1
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="SAR1A")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="SAR1A"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg07661904","cg14213021")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
SAR1A_meth <- MA_mean
colnames(SAR1A_meth) <- "SAR1A"
group_SAR1A_meth <- kmeans(SAR1A_meth,2,nstart = 30)
group_SAR1Ameth_surv <- group_SAR1A_meth$cluster
fit_SAR1Ameth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_SAR1Ameth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_SAR1Ameth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_SAR1A_meth=p_value
p_SAR1A_meth <- signif(p_value,3)


tiff(filename = "SAR1A.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_SAR1Ameth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_SAR1A_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#2
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="NFE2L3")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="NFE2L3"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg13118545","cg14534464","cg16882373","cg18844118")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
NFE2L3_meth <- MA_mean
colnames(NFE2L3_meth) <- "NFE2L3"
group_NFE2L3_meth <- kmeans(NFE2L3_meth,2,nstart = 30)
group_NFE2L3meth_surv <- group_NFE2L3_meth$cluster
fit_NFE2L3meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_NFE2L3meth_surv),
                               data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_NFE2L3meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_NFE2L3_meth=p_value
p_NFE2L3_meth <- signif(p_value,3)


tiff(filename = "NFE2L3.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_NFE2L3meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_NFE2L3_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#3
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="EVI2A")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="EVI2A"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg20368567","cg23352695")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
EVI2A_meth <- MA_mean
colnames(EVI2A_meth) <- "EVI2A"
group_EVI2A_meth <- kmeans(EVI2A_meth,2,nstart = 30)
group_EVI2Ameth_surv <- group_EVI2A_meth$cluster
fit_EVI2Ameth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_EVI2Ameth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_EVI2Ameth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_EVI2A_meth=p_value
p_EVI2A_meth <- signif(p_value,3)


tiff(filename = "EVI2A.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_EVI2Ameth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_EVI2A_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#4
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="SH2B2")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="SH2B2"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg05636499","cg16476940","cg24544047")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
SH2B2_meth <- MA_mean
colnames(SH2B2_meth) <- "SH2B2"
group_SH2B2_meth <- kmeans(SH2B2_meth,2,nstart = 30)
group_SH2B2meth_surv <- group_SH2B2_meth$cluster
fit_SH2B2meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_SH2B2meth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_SH2B2meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_SH2B2_meth=p_value
p_SH2B2_meth <- signif(p_value,3)


tiff(filename = "SH2B2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_SH2B2meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_SH2B2_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#5
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="RP11-440G5.2")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="RP11-440G5.2"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg14211915")]
MA_mean <- as.data.frame(MA_value)
RP11_440G5.2_meth <- MA_mean
colnames(RP11_440G5.2_meth) <- "RP11-440G5.2"
group_RP11_440G5.2_meth <- kmeans(RP11_440G5.2_meth,2,nstart = 30)
group_RP11_440G5.2meth_surv <- group_RP11_440G5.2_meth$cluster
fit_RP11_440G5.2meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_RP11_440G5.2meth_surv),
                                     data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_RP11_440G5.2meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_RP11_440G5.2_meth=p_value
p_RP11_440G5.2_meth <- signif(p_value,3)


tiff(filename = "RP11_440G5.2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片

windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_RP11_440G5.2meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_RP11_440G5.2_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#6
which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="FCGRT")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="FCGRT"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg00586700","cg09293816","cg11847636","cg15528736","cg18081863")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
FCGRT_meth <- MA_mean
colnames(FCGRT_meth) <- "FCGRT"
group_FCGRT_meth <- kmeans(FCGRT_meth,2,nstart = 30)
group_FCGRTmeth_surv <- group_FCGRT_meth$cluster
fit_FCGRTmeth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_FCGRTmeth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_FCGRTmeth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_FCGRT_meth=p_value
p_FCGRT_meth <- signif(p_value,3)


tiff(filename = "FCGRT.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggsurvplot(fit_FCGRTmeth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_FCGRT_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#7

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="PEBP4")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="PEBP4"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg00901161","cg01966160","cg02382109","cg04447123","cg16575322")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
PEBP4_meth <- MA_mean
colnames(PEBP4_meth) <- "PEBP4"
group_PEBP4_meth <- kmeans(PEBP4_meth,2,nstart = 30)
group_PEBP4meth_surv <- group_PEBP4_meth$cluster
fit_PEBP4meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_PEBP4meth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_PEBP4meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_PEBP4_meth=p_value
p_PEBP4_meth <- signif(p_value,3)


tiff(filename = "PEBP4.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_PEBP4meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_PEBP4_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#8

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="EP400")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="EP400"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg05629774","cg08299299","cg12488585","cg23793263")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
EP400_meth <- MA_mean
colnames(EP400_meth) <- "EP400"
group_EP400_meth <- kmeans(EP400_meth,2,nstart = 30)
group_EP400meth_surv <- group_EP400_meth$cluster
fit_EP400meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_EP400meth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_EP400meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_EP400_meth=p_value
p_EP400_meth <- signif(p_value,3)


tiff(filename = "EP400.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_EP400meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_EP400_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#9

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="LINC00982")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="LINC00982"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg00945293","cg01908338","cg04241652","cg05346286","cg07348009","cg07363855","cg08537652","cg08541612","cg11229543","cg12050271","cg15390554","cg25152368","cg26938364")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
LINC00982_meth <- MA_mean
colnames(LINC00982_meth) <- "LINC00982"
group_LINC00982_meth <- kmeans(LINC00982_meth,2,nstart = 30)
group_LINC00982meth_surv <- group_LINC00982_meth$cluster
fit_LINC00982meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_LINC00982meth_surv),
                                  data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_LINC00982meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_LINC00982_meth=p_value
p_LINC00982_meth <- signif(p_value,3)


tiff(filename = "LINC00982.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_LINC00982meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_LINC00982_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#10

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="OSM")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="OSM"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg00017461","cg02026204","cg02597698","cg10467217","cg24727122","cg25666403")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
OSM_meth <- MA_mean
colnames(OSM_meth) <- "OSM"
group_OSM_meth <- kmeans(OSM_meth,2,nstart = 30)
group_OSMmeth_surv <- group_OSM_meth$cluster
fit_OSMmeth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_OSMmeth_surv),
                            data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_OSMmeth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_OSM_meth=p_value
p_OSM_meth <- signif(p_value,3)


tiff(filename = "OSM.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_OSMmeth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_OSM_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#11

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="ECEL1P2")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="ECEL1P2"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg02812891","cg09490371","cg13138089","cg24767131","cg25560398")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
ECEL1P2_meth <- MA_mean
colnames(ECEL1P2_meth) <- "ECEL1P2"
group_ECEL1P2_meth <- kmeans(ECEL1P2_meth,2,nstart = 30)
group_ECEL1P2meth_surv <- group_ECEL1P2_meth$cluster
fit_ECEL1P2meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_ECEL1P2meth_surv),
                                data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_ECEL1P2meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_ECEL1P2_meth=p_value
p_ECEL1P2_meth <- signif(p_value,3)


tiff(filename = "ECEL1P2.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_ECEL1P2meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_ECEL1P2_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#12

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="MIR5687")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="MIR5687"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg26544665")]
#MA_mean <- as.data.frame(apply(MA_value,1,mean))
MIR5687_meth <- as.data.frame(MA_value)
colnames(MIR5687_meth) <- "MIR5687"
group_MIR5687_meth <- kmeans(MIR5687_meth,2,nstart = 30)
group_MIR5687meth_surv <- group_MIR5687_meth$cluster
fit_MIR5687meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_MIR5687meth_surv),
                                data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_MIR5687meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_MIR5687_meth=p_value
p_MIR5687_meth <- signif(p_value,3)


tiff(filename = "MIR5687.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_MIR5687meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_MIR5687_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#13

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="NFU1")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="NFU1"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg03375371","cg10354495","cg15844365","cg20547295")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
NFU1_meth <- MA_mean
colnames(NFU1_meth) <- "NFU1"
group_NFU1_meth <- kmeans(NFU1_meth,2,nstart = 30)
group_NFU1meth_surv <- group_NFU1_meth$cluster
fit_NFU1meth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_NFU1meth_surv),
                             data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_NFU1meth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_NFU1_meth=p_value
p_NFU1_meth <- signif(p_value,3)


tiff(filename = "NFU1.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_NFU1meth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_NFU1_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#14

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="IL1B")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="IL1B"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg07935264","cg18773937","cg23149881")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
IL1B_meth <- MA_mean
colnames(IL1B_meth) <- "IL1B"
group_IL1B_meth <- kmeans(IL1B_meth,2,nstart = 30)
group_IL1Bmeth_surv <- group_IL1B_meth$cluster
fit_IL1Bmeth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_IL1Bmeth_surv),
                             data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_IL1Bmeth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_IL1B_meth=p_value
p_IL1B_meth <- signif(p_value,3)


tiff(filename = "IL1B.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_IL1Bmeth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_IL1B_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()
#15

which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="BGLAP")
MA=methy_tumor_KIRP30knn[which(methy_tumor_KIRP30knn[,"Gene_Symbol"]=="BGLAP"),]
#取平均
MA_value <- methy_KIRC_match[,c("cg11519176","cg12036877","cg17419597","cg23729443","cg24843609")]
MA_mean <- as.data.frame(apply(MA_value,1,mean))
BGLAP_meth <- MA_mean
colnames(BGLAP_meth) <- "BGLAP"
group_BGLAP_meth <- kmeans(BGLAP_meth,2,nstart = 30)
group_BGLAPmeth_surv <- group_BGLAP_meth$cluster
fit_BGLAPmeth_surv <- survfit(Surv(overall_survival,as.numeric(vital_status))~(group_BGLAPmeth_surv),
                              data= KIRC_clinical)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival, as.numeric(vital_status)) ~ group_BGLAPmeth_surv, 
             data = KIRC_clinical)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p_BGLAP_meth=p_value
p_BGLAP_meth <- signif(p_value,3)


tiff(filename = "BGLAP.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggsurvplot(fit_BGLAPmeth_surv, linetype = "strata",
           size=1.6,
           # 不同 strata 使用不同的线性
           legend.labs = c("Cluster1", "Cluster2"), # 设置图例
           palette = c("#ff7473", "#4bc2c5"),  #设置颜色
           ggtheme = theme_bw(), # 设置主题
           font.x = c(15, "black"),
           font.y = c(15, "black"),
           font.legend=c(13,'plain','black'), #设置标签字体
           xlab='Survival time(month)',
           ylab='Survival probability',
           font.tickslab = c(15, "plain", "black"), #设置坐标轴刻度的字体大小和颜色
           # 添加 pvalue
           pval = p_BGLAP_meth, # 显示 pvalue
           pval.size = 6
)
dev.off()

#热图
heat_methy1_15=cbind(SAR1A_meth,NFE2L3_meth,EVI2A_meth,SH2B2_meth,RP11_440G5.2_meth,FCGRT_meth,PEBP4_meth,EP400_meth,LINC00982_meth,OSM_meth,ECEL1P2_meth,MIR5687_meth,NFU1_meth,IL1B_meth,BGLAP_meth)
heat_methy_15=t(heat_methy1_15)[,colnames(miRNA_differentgene_pheatmap)]

methy_differentgene15_pheatmap<-heat_methy_15

tiff(filename = "methy-RRCC.tiff",res = 300,width = 3.6,height = 1.3,units = 'in',compression = c("lzw"))
library(pheatmap)
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
         fontsize = 6.2,
         legend = TRUE,
         annotation_legend = TRUE #不显示图例
)

dev.off()
##############################################
#
# part 8 miRNA target genes                  #
#
###############################################
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/miRNA")
miRNA_fortarget=c('hsa-miR-139-5p',
                  'hsa-miR-21',
                  'hsa-miR-204-5p',
                  'hsa-miR-10b',
                  'hsa-miR-1468',
                  'hsa-miR-130b', 
                  'hsa-miR-625-3p', 
                  'hsa-miR-221-3p', 
                  'hsa-miR-193b',
                  'hsa-miR-9-1',
                  'hsa-miR-9-2',
                  'hsa-miR-9-3',
                  'hsa-miR-671-5p',
                  'hsa-miR-27b',
                  'hsa-miR-365a')
#miRNA_fortarget=c('hsa-miR-4326')
miRNA_fortarget=as.data.frame(miRNA_fortarget)
colnames(miRNA_fortarget)='miRNA'
load("F:/课题代码/human_miRNA_target.Rdata")
miRNA_BRCA_target=merge(miRNA_fortarget, human_miRNA_target, by = 'miRNA', all= F,sort=F)
Target_Gene=unique(miRNA_BRCA_target$`Target Gene`)
setwd("F:\\ecpcs结果\\KIRC下游\\differential gene/miRNA")
write.csv(Target_Gene,file='Target_Gene.csv')





#go富集
setwd("F:\\ecpcs结果\\KIRC下游/go")

enrich_gene=unique(c(Target_Gene,sig_gene_methy_sort_15,sig_gene_mRNA_sort_15[,1]))
write.csv(enrich_gene,"enrich_gene.csv")#前15+15+靶基因

library(ggplot2)
kobas_GO=read.csv("F:\\ecpcs结果\\KIRC下游\\go\\go10.csv")
tiff(filename = "GO.tiff",res = 300,width = 9.5,height = 4,units = 'in',compression = c("lzw")) #输出dpi=300的图片
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


kobas_kegg=read.csv("F:\\ecpcs结果\\KIRC下游\\go\\kegg10.csv")
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


########################################################################################
#
#pathway activity analysis   ECPCS #分两型
#
########################################################################################
library(progeny)
#library(ggbiplot)
library(RColorBrewer)
library(igraph)


#########################################################
# step1 calculate pathway activities for mRNA_KIRP_match #
#########################################################
mRNA_KIRP_match <- as.matrix(t(mRNA_LGG_delmeancv)) #progeny()函数要求输入数据必须为矩阵,行为基因，列为样本
score_pathways <- progeny(mRNA_KIRP_match, scale=TRUE, organism="Human", top = 100, perm = 1)


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


save(differ_pathway_fdr,file='differ_pathway_fdr.Rdata')
write.csv(differ_pathway_fdr,file='differ_pathway_fdr_ecpcs.csv')

#########################################################
# step3 boxplot for average pathway activities in different clusters #
#########################################################

library(ggplot2)
library(ggpubr)

pathway_box <- cbind(score_pathways, group_gaodi)
pathway_box<-as.numeric(pathway_box)
## method 1 for boxplot
mainname=row.names(differ_pathway_fdr)[11] #可以修改g，使箱式图分别显示不同通路结果
mainname
#group_gaodi=as.data.frame(group_gaodi)
TNFα <- as.data.frame(pathway_box[,c("TNFa","group_gaodi")])
colnames(TNFα) <- c("Pathway activity","Group")
VEGF <- as.data.frame(pathway_box[,c("VEGF","group_gaodi")])
colnames(VEGF) <- c("Pathway activity","Group")
Androgen <- as.data.frame(pathway_box[,c("Androgen","group_gaodi")])
colnames(Androgen) <- c("Pathway activity","Group")
EGFR <- as.data.frame(pathway_box[,c("EGFR","group_gaodi")])
colnames(EGFR) <- c("Pathway activity","Group")
NFkβ <- as.data.frame(pathway_box[,c("NFkB","group_gaodi")])
colnames(NFkβ) <- c("Pathway activity","Group")
PI3K <- as.data.frame(pathway_box[,c("PI3K","group_gaodi")])
colnames(PI3K) <- c("Pathway activity","Group")





setwd("F:\\ecpcs结果\\KIRC下游/pathway activity")
setwd("F:\\ecpcs结果\\KIRC下游")
mysize = 6 #p值大小设置
tiff(filename = "TNFα.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(TNFα,x='Group',y='Pathway activity',
          color="black",palette = c("#D25565","#79BEDB"),
          title = "TNFα",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=3,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "myFont",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "myFont",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "VEGF.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(VEGF,x='Group',y='Pathway activity',
          color="black",palette = c("#D25565","#79BEDB"),
          title = "VEGF",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.6,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "Androgen.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(Androgen,x='Group',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "Androgen",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.6,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times New Roman",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times New Roman",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
tiff(filename = "EGFR.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
ggboxplot(EGFR,x='Group',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "EGFR",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=4,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()


tiff(filename = "NFkβ.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(NFkβ,x='Group',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "NFkβ",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=2.6,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()

tiff(filename = "PI3K.tiff",res = 300,width = 4,height = 5,units = 'in',compression = c("lzw")) #输出dpi=300的图片
windowsFonts(myFont = windowsFont("Times New Roman")) ##字体
par(family="myFont") #设置后图里面字体均为Times New Roman
ggboxplot(PI3K,x='Group',y='Pathway activity',
          color="black",fill = c("#D25565","#79BEDB"),
          title = "PI3K",font.title=18,
          font.x = c(18, "black"),font.y = c(18, "black"))+
  stat_compare_means(aes(label=paste0("p=",..p.format..)),label.x=0.7,label.y=6.3,size = mysize)+ #改变label.x/y的数值调整P值在图中的位置
  theme(plot.title=element_text(hjust = 0.5))+ #将图标题居中
  theme(axis.text.x = element_text(colour = "black",family = "Times",size = 14), #设置x轴刻度字体
        axis.text.y = element_text(family = "Times",size = 15)) #设置x,y轴刻度字体,字体样式为“plain"
dev.off()






