addpath(genpath('Diff'));%��Ӽ����غ��������ļ���·��
[v w]=xlsread("pan-cancer.xlsx");
cancer_list=w(:,1);
mRNA_normal_data_all=[];
for i=1:length(v)
    cancer=cancer_list{i};
    load(['./data/',cancer,'/',cancer,'_gene_expression_RNAseq.mat'])
    mRNA_normal_data_all = [mRNA_normal_data_all mRNA_normal_data];
end
for i=1:length(v)
    cancer=cancer_list{i};
    diff_transcription2(cancer,mRNA_normal_data_all);
end