addpath(genpath('co_net'));%��Ӽ����غ��������ļ���·��
[v w]=xlsread("pan-cancer.xlsx");
cancer_list=w(:,1);
for i=1:length(v)
cancer=cancer_list{i};

demo_Meco(cancer);
end