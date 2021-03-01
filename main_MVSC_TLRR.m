%% before running the code, you need download the LinADMM library from:https://github.com/canyilu/LibADMM

%% load data
clear
clc
load('EYaleB10_mtv.mat') % the input multiple view dataw 
load('YaleB_IJCV') % the initial input muptiple view similarity matrix. This similarity matrices was generated from 
% the work <On Unifying Multi-View Self-Representations for Clustering by Tensor Multi-Rank Minimization, 
% International Journal of Computer Vision, 2018>
%% 
addpath(genpath('..\LibADMM-master'))
addpath('.\misc')
addpath('.\ClusteringMeasure')

%% parameter settings 
para.lambda=15;
para.w1=0.4;  
para.alpha=9; 

% 
% The input similarity tensor construciton 
Wtensor(:,:,1)=abs(Z{1})+abs(Z{1}');
Wtensor(:,:,2)=abs(Z{2})+abs(Z{2}');
Wtensor(:,:,3)=abs(Z{3})+abs(Z{3}');


[L,E]=fun_MVSC_TLRR(Wtensor,para);




%%  construct the output similarity tensor
close all
S=abs(L(:,:,1))+abs(L(:,:,2))+abs(L(:,:,3));

S=S-diag(diag(S));
imagesc(S)
colormap(jet)
truesize([300 300]);



%% compute ACC and NMI
cls_num = length(unique(gt));
for ii=1:20
C = SpectralClustering(S,cls_num);
[A nmi(ii) avgent] = compute_nmi(gt,C);
ACC(ii) = Accuracy(C,double(gt));
end
disp('ACC,NMI')
disp([mean(ACC),mean(nmi)])