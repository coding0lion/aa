addpath(genpath('co_net'));%添加计算熵函数所在文件夹路径
[v w]=xlsread("pan-cancer.xlsx");
cancer_list=w(:,1);
for i=1:length(v)
cancer=cancer_list{i};

demo_Meco(cancer);
end