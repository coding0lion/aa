function [gene,data]=loda_data_Meco(cancer)
load(['./data/',cancer,'/',cancer,'_met.mat'])%load����������
gene=met_gene;
data=met_patient_data;
end