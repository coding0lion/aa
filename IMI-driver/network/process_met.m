[v w]=xlsread("pan-cancer.xlsx");
cancer_list=w(:,1);
for i=1:length(v)
cancer=cancer_list{i};
[data word]=xlsread(['./data/raw/',cancer,'_Me.xlsx']);
met_gene=word(2:end,1);
sample=word(1,2:end);
ID=sample;
id_num=[];
%得到序列的后两位
for i=1:length(ID)
id_num(i,1)=str2double(ID{1,i}(end-1:end));
end
%正常样本
%索引
normal_index=[];
met_normal_ID=[];
met_normal_data=[];
normal_index=find(id_num>=10);

met_normal_ID=ID(normal_index);
met_normal_data=data(:,normal_index);

patient_index=[];
met_patient_ID=[];
met_patient_data=[];
%病人样本id
patient_index=find(id_num<10);

met_patient_ID=ID(patient_index);
met_patient_data=data(:,patient_index);

save(['./output/',cancer,'/',cancer,'_met.mat'], 'met_gene', 'met_patient_data', 'met_patient_ID','met_normal_data','met_normal_ID','-v7.3')
end