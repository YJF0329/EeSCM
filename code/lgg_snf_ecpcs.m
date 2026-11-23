clc
clear
close all
%quant10
load('methyLGGdelmeancv.mat')
load('miRNALGGdelmeancv.mat')
load('mRNALGGdelmeancv.mat')
  
  X_mRNA=mRNALGGdelmeancv;
  X_miRNA=miRNALGGdelmeancv;
  X_methy=methyLGGdelmeancv;
  
  
  n=size(X_mRNA);%维度
  m=n(1,1);%m样本数
%%snf
  addpath(genpath('./SNFmatlab'));
  SXmRNA=Standard_Normalization(X_mRNA);%行是样本，列是基因
  SXmiRNA=Standard_Normalization(X_miRNA);
  SXmethy=Standard_Normalization(X_methy);
  K = 30;		% number of neighbors, usually (10~30)
  alpha = 0.5; % hyperparameter, usually (0.3~0.8)
  T= 20;
  Dist1_snf = dist2(SXmRNA,SXmRNA);
  Dist2_snf = dist2(SXmiRNA,SXmiRNA);
  Dist3_snf = dist2(SXmethy,SXmethy);
 W1_snf = affinityMatrix(Dist1_snf, K, alpha);
 W2_snf = affinityMatrix(Dist2_snf, K, alpha);
 W3_snf = affinityMatrix(Dist3_snf, K, alpha);
 W_snf = SNF({W1_snf,W2_snf,W3_snf}, K, T);
 group2_snf = SpectralClustering(W_snf,2);
 group3_snf = SpectralClustering(W_snf,3);
  group4_snf = SpectralClustering(W_snf,4);
  group5_snf = SpectralClustering(W_snf,5);
group6_snf = SpectralClustering(W_snf,6);
group7_snf = SpectralClustering(W_snf,7);
group8_snf = SpectralClustering(W_snf,8);
group9_snf = SpectralClustering(W_snf,9);
group10_snf = SpectralClustering(W_snf,10);
[K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W_snf, [2:5]);
fprintf('The best number of clusters according to eigengap is %d\n', K1);
fprintf('The best number of clusters according to rotation cost is %d\n', K2);
%%ne-snf
%永清师姐在融合矩阵W降噪
 addpath(genpath('./ne降噪Network_Enhancement-master'));
 W_NE=Network_Enhancement(W_snf);
 group2_NE = SpectralClustering(W_NE,2);
group3_NE = SpectralClustering(W_NE,3);
group4_NE = SpectralClustering(W_NE,4);
group5_NE = SpectralClustering(W_NE,5);
 group6_NE = SpectralClustering(W_NE,6);
 group7_NE = SpectralClustering(W_NE,7);
 group8_NE = SpectralClustering(W_NE,8);
 group9_NE = SpectralClustering(W_NE,9);
 group10_NE = SpectralClustering(W_NE,10);

 %all(all(W1_skf==W1_skf==1)==1)
 %%skf
 %用skf的标准化
 addpath(genpath('./SKFCancer-subtypes-master'));
 SXmRNA1 = zscore(X_mRNA.');%转置为 列是样本
SXmiRNA2 = zscore(X_miRNA.');
SXmethy3 = zscore(X_methy.');
W1_skf = kernel_eucl(SXmRNA1,2);  
W2_skf = kernel_eucl(SXmiRNA2,2);  
W3_skf = kernel_eucl(SXmethy3,2); 
K = 30;%number of neighbors, usually (10~30)
alpha = 0.5; %hyperparameter, usually (0.3~0.8)
T = 20; %Num
W_skf = SKF({W1_skf,W2_skf,W3_skf}, K, T,alpha);
group2_skf = SpectralClustering(W_skf,2);
group3_skf = SpectralClustering(W_skf,3);
group4_skf = SpectralClustering(W_skf,4);
group5_skf = SpectralClustering(W_skf,5);  
group6_skf = SpectralClustering(W_skf,6);
group7_skf = SpectralClustering(W_skf,7);
group8_skf = SpectralClustering(W_skf,8);
group9_skf = SpectralClustering(W_skf,9);
group10_skf = SpectralClustering(W_skf,10);
save group2_skf.mat group2_skf;
save group3_skf.mat group3_skf;
save group4_skf.mat group4_skf;
save group5_skf.mat group5_skf;
save group6_skf.mat group6_skf;
save group7_skf.mat group7_skf;
save group8_skf.mat group8_skf;
save group9_skf.mat group9_skf;
save group10_skf.mat group10_skf;

%csvwrite('W_skf2.csv',W_skf)
%csvwrite('group2_skf2.csv',group2_skf)
%csvwrite('group3_skf2.csv',group3_skf)
%csvwrite('group4_skf2.csv',group4_skf)
%csvwrite('group5_skf2.csv',group5_skf)
 %%hopes
 addpath(genpath('./HOPES-master'));
K = 30;%neighbor number
 W1_hopes = affinityMatrix(Dist1_snf, K );
 W2_hopes = affinityMatrix(Dist2_snf, K );
 W3_hopes = affinityMatrix(Dist3_snf, K);
 beta = 0.25;
 alpha = 0.5;
 W_hopes = HOPES({W1_hopes,W2_hopes,W3_hopes},K,alpha,beta);
 group2_hopes = SpectralClustering(W_hopes,2);
group3_hopes = SpectralClustering(W_hopes,3);
group4_hopes = SpectralClustering(W_hopes,4); 
group5_hopes = SpectralClustering(W_hopes,5);
 group6_hopes = SpectralClustering(W_hopes,6);
 group7_hopes = SpectralClustering(W_hopes,7);
 group8_hopes = SpectralClustering(W_hopes,8);
 group9_hopes = SpectralClustering(W_hopes,9);
 group10_hopes = SpectralClustering(W_hopes,10);

%%Movics
addpath(genpath('./RunRcode'));
save W_snf.mat W_snf;
save group2_snf.mat group2_snf;
save group3_snf.mat group3_snf;
save group4_snf.mat group4_snf;
save group5_snf.mat group5_snf;
save group6_snf.mat group6_snf;
save group7_snf.mat group7_snf;
save group8_snf.mat group8_snf;
save group9_snf.mat group9_snf;
save group10_snf.mat group10_snf;
save W_NE.mat W_NE;
save group2_NE.mat group2_NE;
save group3_NE.mat group3_NE;
save group4_NE.mat group4_NE;
save group5_NE.mat group5_NE;
save group6_NE.mat group6_NE;
save group7_NE.mat group7_NE;
save group8_NE.mat group8_NE;
save group9_NE.mat group9_NE;
save group10_NE.mat group10_NE;
save W_skf.mat W_skf;
save group2_skf.mat group2_skf;
save group3_skf.mat group3_skf;
save group4_skf.mat group4_skf;
save group5_skf.mat group5_skf;
save group6_skf.mat group6_skf;
save group7_skf.mat group7_skf;
save group8_skf.mat group8_skf;
save group9_skf.mat group9_skf;
save group10_skf.mat group10_skf;

save W_hopes.mat W_hopes;
save group2_hopes.mat group2_hopes;
save group3_hopes.mat group3_hopes;
save group4_hopes.mat group4_hopes;
save group5_hopes.mat group5_hopes;
save group6_hopes.mat group6_hopes;
save group7_hopes.mat group7_hopes;
save group8_hopes.mat group8_hopes;
save group9_hopes.mat group9_hopes;
save group10_hopes.mat group10_hopes;


%%ECPCS
addpath(genpath('./ECPCS'));
 baseCLs=zeros(m,4);%n要换
 baseCLs(:,1)=group5_snf;
 baseCLs(:,2)=group4_NE;
 baseCLs(:,3)=group4_skf;
 baseCLs(:,4)=group2_hopes;
 [bcs, baseClsSegs] = getAllSegs(baseCLs);%bcs:将类总排序；baseClsSegs：类数*样本数 在这一类中有这一样本为1，否则0.
clsSim = full(simxjac(baseClsSegs)); % Build the cluster similarity matrix by Jaccard coefficient.  类数*类数
t=20;
clsSimRW = computePTS_II(clsSim, t); % Perform random walks and obtain a new cluster-wise similarity matrix.类簇数*类簇数
ECA = getECA(bcs,clsSimRW); %样本数*样本数
 group2_ECPCS_HC = ECPCS_HC(baseCLs, 2);
 save ECA.mat ECA;

 group3_ECPCS_HC = ECPCS_HC(baseCLs, 3);
 group4_ECPCS_HC = ECPCS_HC(baseCLs, 4);
 %csvwrite('group4_ECPCS_HC1.csv',group4_ECPCS_HC)
 group5_ECPCS_HC = ECPCS_HC(baseCLs, 5);
 group6_ECPCS_HC = ECPCS_HC(baseCLs, 6); 
 group7_ECPCS_HC = ECPCS_HC(baseCLs, 7);
 group8_ECPCS_HC = ECPCS_HC(baseCLs, 8);
 group9_ECPCS_HC = ECPCS_HC(baseCLs, 9);
 group10_ECPCS_HC = ECPCS_HC(baseCLs, 10);
save group2_ECPCS_HC.mat group2_ECPCS_HC;
save group3_ECPCS_HC.mat group3_ECPCS_HC;
save group4_ECPCS_HC.mat group4_ECPCS_HC;
save group5_ECPCS_HC.mat group5_ECPCS_HC;
save group6_ECPCS_HC.mat group6_ECPCS_HC;
save group7_ECPCS_HC.mat group7_ECPCS_HC;
save group8_ECPCS_HC.mat group8_ECPCS_HC;
save group9_ECPCS_HC.mat group9_ECPCS_HC;
save group10_ECPCS_HC.mat group10_ECPCS_HC;
