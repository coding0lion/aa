[~,cancer_t] = xlsread('cancer.xlsx');
for i = 1:length(cancer_t) 
cancer = cancer_t{i};
clear word data
load([ './output/',cancer ,'/' ,cancer,'_Dif_transcription.mat'])
Dif_transcription(isnan(Dif_transcription))=0;
[data word]=xlsread(['./data/',cancer,'/' ,cancer,'_Dif_Methylation.xlsx']);
m_gene = word(2:end,1);
m_d = zeros(length(gene),1);
for i = 1:length(m_gene)
    index =find(strcmp(m_gene{i},gene));
    m_d(index) = data(i);
end
dif = [m_d Dif_transcription];
save_path = ['./output/dif_data/Dif_',cancer,'.csv'];
writetable(table(dif),save_path);
end
