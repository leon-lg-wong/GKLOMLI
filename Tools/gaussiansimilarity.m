function [DS,RS]=gaussiansimilarity(A)

[disease_num,rna_num]=size(A);


Gaussian_Dis = zeros(disease_num,disease_num);
pare_a=0;
sum=0;
temp=0;

for i=1:disease_num
    temp=norm(A(i,:));
    sum=sum+temp^2;
end
pare_a=1/(sum/disease_num);

for i=1:disease_num
    
    for j=1:disease_num
        Gaussian_Dis(i,j)=exp(-pare_a*(norm(A(i,:)-A(j,:))^2));
    end
end
DS=Gaussian_Dis;

Gaussian_miR = zeros(rna_num,rna_num);
pare_b=0;
sum=0;
temp=0;

for i=1:rna_num
    temp=norm(A(:,i));
    sum=sum+temp^2;
end
pare_b=1/(sum/rna_num);

for i=1:rna_num
    
    for j=1:rna_num
        Gaussian_miR(i,j)=exp(-pare_b*(norm(A(:,i)-A(:,j))^2));
    end
end
RS=Gaussian_miR;
