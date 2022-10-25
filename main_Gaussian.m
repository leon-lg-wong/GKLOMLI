clear
tic

load /Dataset/miRNA-lnRNA_interaction/interaction.mat
addpath("./Tools/")

num_Known_Association=length(find(interaction==1));
n_test=100;
k_fold=5;
overallauc=zeros(n_test,1);
Random_order=zeros(num_Known_Association,n_test);
for cv=1:n_test
    Random_order(:,cv)=crossvalind("KFOLD",num_Known_Association,k_fold);
end

parfor cv=1:n_test
    if overallauc(cv,1)~=0
        continue;
    end
    [POSITION(cv,:)]=Method_Gaussian(interaction, k_fold,Random_order(:,cv),0.018);
	
	[fpr(cv,:),tpr(cv,:),overallauc(cv,1)] = positiontooverallauc(POSITION(cv,:),interaction,k_fold);	
end

Auc_avg=mean(overallauc(:,1));
AUC_std=std(overallauc(:,1));
if exist('Results','dir')==0
    mkdir("Results");
end
save Results/Gaussian_results.mat POSITION Auc_avg AUC_std Random_order fpr tpr
toc