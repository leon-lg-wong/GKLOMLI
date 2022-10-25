function [position]=Method_Single(Network_AB, k_fold,index_kf,alpha_A)
ori_interaction=Network_AB;
index_ori=find(ori_interaction==1);

T=1;
for ccv = 1 : k_fold
    index_test=index_ori(index_kf==ccv);    
    num_test=length(index_test);
    Target_network=ori_interaction;
    Target_network(index_test)=0;
    
    train=Target_network;
    
    temp=train; % Integration matrix \mathbf{A}
    
    loMatrixA = (1/alpha_A*eye(size(temp'*temp)) + ...
        temp' * temp)\temp' * temp; % Computing the weighting matrix
    loMatrixA = temp * loMatrixA;
    F = loMatrixA;
    
    finalscore=F(index_test);
    % make the score of seed miRNA-lncRNAinteractions as zero
    
    F(index_ori(index_kf~=ccv))=-99999999999999999999;
    for qq=1:num_test
          ll1=size(find(F>=finalscore(qq)),1);
        ll2=size(find(F>finalscore(qq)),1);
        position(1,T)=ll2+1+(ll1-ll2-1)/2;
        T=T+1;
    end
    
end
