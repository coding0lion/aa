[v w]=xlsread("diff_cancer_29.xlsx");
cancer_list=w(:,1);
for i=2:length(w)
cancer=cancer_list{i};
addpath(genpath('pagerank'));%��Ӽ����غ��������ļ���·��
%all_pagerank(cancer);
all_pagerank_nodiff(cancer);
end
