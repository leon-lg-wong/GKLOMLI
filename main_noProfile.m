
tic
load ./Dataset/miRNA-lnRNA_interaction/interaction.mat
addpath("./Tools/")

num_Known_Association=length(find(interaction==1));
n_test=100;
k_fold=5;
overallauc=zeros(n_test,1);
Random_order=zeros(num_Known_Association,n_test);
for cv=1:n_test
    Random_order(:,cv)=crossvalind("KFOLD",num_Known_Association,k_fold);
end
alph=0.001:0.0005:0.02;
alph=alph';
alph_list=zeros(n_test,1);
for cv=1:n_test
    if overallauc(cv,1)~=0
        continue;
    end
    
    tmp_auc=zeros(length(alph),1);
    
    parfor alph_i = 1: length(alph)
        
        [tmp_POSITION(alph_i,:)]=Method_Single(interaction, k_fold,Random_order(:,cv),alph( alph_i ));
        [tmp_fpr(alph_i,:),tmp_tpr(alph_i,:),tmp_auc(alph_i)] = positiontooverallauc(tmp_POSITION(alph_i,:),interaction,k_fold);
    end
    index_tmp_auc=find(tmp_auc == max(tmp_auc));
    fpr(cv,:)=tmp_fpr(index_tmp_auc(1),:);
    tpr(cv,:)=tmp_tpr(index_tmp_auc(1),:);
    overallauc(cv,1) = tmp_auc(index_tmp_auc(1));
    POSITION(cv,:) = tmp_POSITION(index_tmp_auc(1),:);
    alph_list(cv)=alph( index_tmp_auc(1) );
end
Auc_avg=mean(overallauc(:,1));
AUC_std=std(overallauc(:,1));
if exist('Results','dir')==0
    mkdir("Results");
end
save ./Results/NoProfile_results.mat POSITION overallauc Auc_avg AUC_std Random_order fpr tpr alph_list
toc
