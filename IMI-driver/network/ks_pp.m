[v w]=xlsread("diff_cancer_29.xlsx");
cancer_list=w(:,1);
ks_p=[];
for i=2:length(cancer_list)
cancer =  cancer_list{i};  
load(['./output/',cancer,'/',cancer,'_rank(cgc).mat'])%load基因表达数据
ks_p = [ks_p;p_f];
end
save('ks_p.mat','ks_p');
xlswrite('ks_p.xlsx',ks_p);