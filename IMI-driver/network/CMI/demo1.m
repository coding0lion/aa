
% load('TCGA_BRAC_mRNA.mat');
[dmfs_e, dmfs_time ,sample_ID ,ID ,data, probeSet]=loda_data();
%PPI��������
load('g_g.mat');

[ma, phenotype]=mapping(dmfs_e,dmfs_time,sample_ID,ID,data);
%���в��Եĸ�ֵ
% ma=randn(30,length(phenotype));
% gene=probeSet(1:30,1);
% for i=1:1000
%     a=randi(50);
%     b=randi(70);
%     candidate_pair{i,1}=probeSet{a,1};
%     candidate_pair{i,2}=probeSet{b,1};
% end

%��ɢ�������л���Ϣ����
% phenotype(1:300,1)=0;
% phenotype(301:593,1)=1;
%%%%phenotype��n*1��n�������ı���(0����1)
%%%%%ʵ������

gene=probeSet;
%����ppi������м���
candidate_pair(:,1)=g_g.X1;
candidate_pair(:,2)=g_g.X2;

permute_type=1;
p_cut=0.05;

tic;
%%�������������������Ҫ����GR_CMI.m
Pair = GR_CMI(ma,gene,phenotype,permute_type,p_cut,candidate_pair);
%Pair = GR_CMI(ma,gene,phenotype,permute_type,p_cut);%�����������ppi���磬��ֱ��ɾȥcandidate_pair
toc;
