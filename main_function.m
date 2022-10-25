clear
tic
load Function_data/lncrna_function_similairty_matrix.mat
load Function_data/mirna_function_similairty_matrix.mat
load Function_data/invalid_lnc_function.mat
load Function_data/invalid_mi_function.mat
load /Dataset/miRNA-lnRNA_interaction/interaction.mat
addpath("./Tools/")

num_Known_Association=length(find(interaction==1));
Sim1=lnc_function_similarity;
Sim2=mi_function_similarity;

%avoid the element of NaN
Sim1(invalid_lnc_function,:) = 0;
Sim1(:,invalid_lnc_function) = 0;
Sim2(invalid_mi_function,:) = 0;
Sim2(:,invalid_mi_function) = 0;

mean_lnc=mean(Sim1(:));
mean_mi=mean(Sim2(:));
Sim1(invalid_lnc_function,:) = mean_lnc;
Sim1(:,invalid_lnc_function) = mean_lnc;
Sim2(invalid_mi_function,:) = mean_mi;
Sim2(:,invalid_mi_function) = mean_mi;

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

    POSITION(cv,:)=Method_SP(interaction,Sim1, Sim2, k_fold,Random_order(:,cv),0.018);

	
	[fpr(cv,:),tpr(cv,:),overallauc(cv,1)] = positiontooverallauc(POSITION(cv,:),interaction,k_fold);	
end

Auc_avg=mean(overallauc(:,1));
AUC_std=std(overallauc(:,1));
if exist('Results','dir')==0
    mkdir("Results");
end
save Results/function_results.mat POSITION Auc_avg AUC_std Random_order fpr tpr
toc