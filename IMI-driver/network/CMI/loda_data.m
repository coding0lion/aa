function [dmfs_e dmfs_time sample_ID ID data, proSet]=loda_data(cancer)

load(['./data/',cancer,'/',cancer,'_clinic.mat']);%�ٴ�����
load(['./data/',cancer,'/',cancer,'_gene_expression_RNAseq.mat'])%load����������
dmfs_e=patient_dmfs_e;
dmfs_time=patient_dmfs_time;
sample_ID=patient_ID;%�ٴ���������id
ID=mRNA_patient_ID;%ת¼������id
data=mRNA_patient_data;
proSet=mRNA_gene;

end
