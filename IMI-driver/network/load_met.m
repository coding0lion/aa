[v w]=xlsread("pan-cancer.xlsx");
cancer_list=w(:,1);
for i=1:length(v)
cancer=cancer_list{i};
load(['./data/raw/',cancer,'_met.mat'])

save(['./data/',cancer,'/',cancer,'_met.mat'], 'met_gene', 'met_patient_data', 'met_patient_ID','met_normal_data','met_normal_ID','-v7.3')
end
