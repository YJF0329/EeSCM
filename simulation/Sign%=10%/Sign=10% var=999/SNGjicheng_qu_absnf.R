library("ClassDiscovery")
library("R.matlab")
library(SNFtool)
#设置为当前代码路径
setwd("F:\\课题代码")
#载入matlab数据
W_snf = readMat("W_snf.mat")#列表
W_snf=W_snf[["W.snf"]]#矩阵
group_snf = readMat("group_snf.mat")
group_snf=as.integer(group_snf[["group.snf"]])

W_NE = readMat("W_NE.mat")
W_NE =W_NE[["W.NE"]]
group_NE = readMat("group_NE.mat")
group_NE =as.integer(group_NE[["group.NE"]])


W_skf = readMat("W_skf.mat")
W_skf =W_skf[["W.skf"]]
group_skf = readMat("group_skf.mat")
group_skf =as.integer(group_skf[["group.skf"]])



W_hopes = readMat("W_hopes.mat")
W_hopes =W_hopes[["W.hopes"]]
group_hopes = readMat("group_hopes.mat")
group_hopes =as.integer(group_hopes[["group.hopes"]])


#
moic1.res.list=NULL
samID=as.character(1:200)
moic1.res.list[["SNF"]][["fit"]]=W_snf
clust.res.snf=data.frame(samID,clust=group_snf)
moic1.res.list[["SNF"]][["clust.res"]]=clust.res.snf
moic1.res.list[["SNF"]][["mo.method"]]="SNF"

moic1.res.list[["NE"]][["fit"]]=W_NE
clust.res.NE=data.frame(samID,clust=group_NE)
moic1.res.list[["NE"]][["clust.res"]]=clust.res.NE
moic1.res.list[["NE"]][["mo.method"]]="NE"


moic1.res.list[["skf"]][["fit"]]=W_skf
clust.res.skf=data.frame(samID,clust=group_skf)
moic1.res.list[["skf"]][["clust.res"]]=clust.res.skf
moic1.res.list[["skf"]][["mo.method"]]="skf"



moic1.res.list[["hopes"]][["fit"]]=W_hopes
clust.res.hopes=data.frame(samID,clust=group_hopes)
moic1.res.list[["hopes"]][["clust.res"]]=clust.res.hopes
moic1.res.list[["hopes"]][["mo.method"]]="hopes"



moic.res.list=moic1.res.list
moic.res <- lapply(moic.res.list, function(x) x$clust.res)
N.clust <- length(unique(moic.res[[1]]$clust))
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
coca.moic <- cutree(hcs, N.clust)#集成后的聚类
writeMat("group_movics.mat",group_movics=coca.moic)



#将相似性矩阵转为距离矩阵再求欧式距离再层次  变为去掉求欧式距离
hcs1 <- hclust(as.dist(as.matrix(1-similarity.matrix)), "average")

coca.moic1 <- cutree(hcs1, N.clust)#集成后的聚类
writeMat("group_movics1.mat",group_movics1=coca.moic1)

#直接对相似矩阵谱聚类
group_movics2 = SNFtool::spectralClustering(as.matrix(similarity.matrix),N.clust)
writeMat("group_movics2.mat",group_movics2=group_movics2)
