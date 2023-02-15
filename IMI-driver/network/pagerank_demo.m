[v w]=xlsread("diff_cancer_29.xlsx");
cancer_list=w(:,1);
for i=2:length(w)
cancer=cancer_list{i};
addpath(genpath('pagerank'));%添加计算熵函数所在文件夹路径
%all_pagerank(cancer);
all_pagerank_nodiff(cancer);
end
