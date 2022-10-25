function [fpr,tpr,overallauc]=positiontooverallauc(position,interaction,kf)

[n,m]=size(interaction);

pp=length(find(interaction==1));
for k=1:m*n-floor(pp/kf)*(kf-1)
    tp=0;
    for t=1:pp
        if position(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    if k<m*n-pp+floor(pp/kf)+1
        fp=k*pp-tp;
    else fp=floor(pp/kf)*(kf-1)*(m*n-pp+floor(pp/kf))+(pp-floor(pp/kf)*(kf-1))*k-tp;
    end
    fpr(1,k)=fp/(floor(pp/kf)*(kf-1)*(m*n-pp+floor(pp/kf)-1)+(pp-floor(pp/kf)*(kf-1))*(m*n-floor(pp/kf)*(kf-1)-1));
end
% plot(fpr,tpr)
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-floor(pp/kf)*(kf-1)
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc=sum(area);
end