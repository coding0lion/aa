%[v w]=xlsread("pan-cancer.xlsx");
%cancer_list=w(:,1);
%for i=1:length(v)
%cancer=cancer_list{i};
%addpath(genpath('CMI'));%���Ӽ����غ��������ļ���·��
%demo_cmi(cancer);

%addpath(genpath('co_net'));%���Ӽ����غ��������ļ���·��
%demo_co(cancer);

% 
% addpath(genpath('synergistic_network'));%���Ӽ����غ��������ļ���·��
% dome_sy(cancer);

cancer='BRCA';
addpath(genpath('sm_net'));
demo_sm(cancer);

%addpath(genpath('ceRNA'));%���Ӽ����غ��������ļ���·��
%demo_ceRNA(cancer);

%addpath(genpath('co_net'));
%demo_Meco(cancer);
%end
