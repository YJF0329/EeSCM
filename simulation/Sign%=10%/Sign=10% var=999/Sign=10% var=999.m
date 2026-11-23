clc
clear
close all
%% 生成数据
n=200;
p=1000;
info_num=p*0.1;
noise_num=p-info_num;
ns=50; % numbers for each subtype

truelabel = kron([1 2 3 4], ones(1, ns));
label=truelabel';%行向量变列向量
C=4;



num_sample=1;

value_var1=9;                                                                                                                                                     
value_var2=9;
value_var3=9;
result_NMI999=zeros(num_sample,16 );
%result_NMIecpcs=zeros(num_sample,16 );%ecpcs中的NMI

tic
for j = 1:num_sample

   rng(j);%设置随机数种子

  X_mRNA=zeros(n,p);
  X_mRNA(1:50,1:info_num)=1;
  X_mRNA(151:200,1:info_num)=3;
  mu=(ones(1,p)*0);%p个0
  Sigma1 =diag((ones(1,p)*value_var1));%对角线为value_var1
  X_mRNA=X_mRNA+mvnrnd(mu,Sigma1,n);
 
  X_miRNA=zeros(n,p);
  X_miRNA(51:100,1:info_num)=2;
  X_miRNA(151:200,1:info_num)=3;
  Sigma2 =diag((ones(1,p)*value_var2));
  X_miRNA=X_miRNA+mvnrnd(mu,Sigma2,n);
  
  X_methy=zeros(n,p);
  X_methy(1:100,1:info_num)=2; 
  X_methy(101:150,1:info_num)=1;
  X_methy(151:200,1:info_num)=3;
  Sigma3 =diag((ones(1,p)*value_var3));
  X_methy=X_methy+mvnrnd(mu,Sigma3,n);
  
  %%snf
  addpath(genpath('./SNFmatlab'));
  SXmRNA=Standard_Normalization(X_mRNA);
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
 group_snf = SpectralClustering(W_snf,C);
 NMI_snf = Cal_NMI(group_snf, label);
 result_NMI999(j,1)=NMI_snf;
 
 
%%ne-snf
%永清师姐在融合矩阵W降噪
 addpath(genpath('./ne降噪Network_Enhancement-master'));
 W_NE=Network_Enhancement(W_snf);
 group_NE = SpectralClustering(W_NE,C);
 NMI_NE = Cal_NMI(group_NE, label);
 result_NMI999(j,2)=NMI_NE;

 
 %%skf
 %用skf的标准化
 addpath(genpath('./SKFCancer-subtypes-master'));
 SXmRNA1 = zscore(X_mRNA.');
SXmiRNA2 = zscore(X_miRNA.');
SXmethy3 = zscore(X_methy.');
W1_skf = kernel_eucl(SXmRNA1,2);  
W2_skf = kernel_eucl(SXmiRNA2,2);  
W3_skf = kernel_eucl(SXmethy3,2); 
K = 30;%number of neighbors, usually (10~30)
alpha = 0.5; %hyperparameter, usually (0.3~0.8)
T = 20; %Num
W_skf = SKF({W1_skf,W2_skf,W3_skf}, K, T,alpha);
group_skf = SpectralClustering(W_skf,C);
 NMI_skf = Cal_NMI(group_skf, label);
 result_NMI999(j,3)=NMI_skf;

  
 %%hopes
 addpath(genpath('./HOPES-master'));
K = 30;%neighbor number
 W1_hopes = affinityMatrix(Dist1_snf, K );
 W2_hopes = affinityMatrix(Dist2_snf, K );
 W3_hopes = affinityMatrix(Dist3_snf, K);
 beta = 0.25;
 alpha = 0.5;
  
 W_hopes = HOPES({W1_hopes,W2_hopes,W3_hopes},K,alpha,beta);
 group_hopes = SpectralClustering(W_hopes,C);
 NMI_hopes = Cal_NMI(group_hopes, label);
 result_NMI999(j,4)=NMI_hopes;

 

%%Movics
addpath(genpath('./RunRcode'));
save W_snf.mat W_snf;
save group_snf.mat group_snf;
save W_NE.mat W_NE;
save group_NE.mat group_NE;

save W_skf.mat W_skf;
save group_skf.mat group_skf;

save W_hopes.mat W_hopes;
save group_hopes.mat group_hopes;

Rpath = 'D:\sougouxiazai\R-4.1.2\bin'; 
RscriptFileName = 'F:\课题代码\SNGjicheng_qu_absnf.R'; 
RunRcode(RscriptFileName, Rpath); 
load('group_movics.mat')
group_movics=cast(group_movics,'double');
NMI_movics= Cal_NMI(group_movics, label);
result_NMI999(j,5)=NMI_movics;



addpath(genpath('./ECPCS'));
 baseCLs=zeros(200,4);
 baseCLs(:,1)=group_snf;
 baseCLs(:,2)=group_NE;
 baseCLs(:,3)=group_skf;
 baseCLs(:,4)=group_hopes;
 
 group_ECPCS_MC = ECPCS_MC(baseCLs, C);
 NMI_ECPCS_MC = NMImax(group_ECPCS_MC,label);
 NMI_ECPCS_MC1 = Cal_NMI(group_ECPCS_MC,label);
 group_ECPCS_HC = ECPCS_HC(baseCLs, C);
 NMI_ECPCS_HC = NMImax(group_ECPCS_HC,label);
 NMI_ECPCS_HC1 = Cal_NMI(group_ECPCS_HC,label);

 result_NMI999(j,7)=NMI_ECPCS_MC1;
 result_NMI999(j,6)=NMI_ECPCS_HC1;

 
 
 
 end
NMI_mean=mean(result_NMI999,1);
csvwrite('result_2.5_999_NMI.csv',result_NMI999)
csvwrite('NMI_mean_25_999.csv',NMI_mean)
NMI_std=std(result_NMI999);
csvwrite('NMI_std_25_999.csv',NMI_std)
toc