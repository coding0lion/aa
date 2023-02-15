%rank_coTRΪת¼����򹲱������rank���У�rank_coMEΪ�׻�������������rank����
% clear;
% clc;
function all_pagerank(cancer)
load(['./data/',cancer,'/',cancer,'_gene_expression_RNAseq.mat'])
% cgc=readtable('tier1_2.txt');
% save('cgc.mat','cgc');
load('cgc.mat')%���а�֢��ͬ������
cgc_gene=cgc.GeneSymbol;
cgc_lab=cgc.Tier;
know_gene=cgc_gene(find(cgc_lab==1));
know_lab=cgc_lab(find(cgc_lab==1));

gene=mRNA_gene;
gene=sortrows(gene,1);
%%
%���ǩ
for i=1:length(gene)
    index=find(strcmp(gene{i},know_gene)==1);
    if index
    gene_labs(i,1)=cgc_lab(index);
    else
    gene_labs(i,1)=-1;
    end
end
%���ñߺͲ������Ȩ��
w=0.6;
e=0.4;
%%
%ת¼��
load(['./output/',cancer,'/','co_net.mat'])
%load('co_net.mat')
d=0.85;
pair_co=pair_co;
[adj_m,weight]=create_adjacency_matrix_weight(gene,pair_co);
adj_m(logical(eye(size(adj_m))))=0; 
ex=ones(length(gene),1);
r=Generank(adj_m,ex,d);
r_co=num2cell(r);
relt_co=[gene r_co];
[h_coTR,p_coTR]=rank_test(relt_co,gene_labs);%����
relt_co=sortrows(relt_co,-2);%���н�������
rank=(1:length(relt_co))';
rank_coTR=[relt_co num2cell(rank)];
rank_coTR=sortrows(rank_coTR,1);%���յ���������
%%
%�׻������򹲱������
load(['./output/',cancer,'/','Meco_net.mat'])
% load('Meco_net.mat')%output
load(['./data/',cancer,'/',cancer,'_met.mat'])
% load('BRCA_met.mat')%�����˺������˷��� data
%[V,W]=xlsread(['./data/',cancer,'/',cancer,'_Dif_Methylation.xlsx']);
% save('BRCA_Dif_methylation','Dif_methylation');
%load('BRCA_Dif_methylation.mat')
%֮ǰ˵���ڽӾ�����ΪmRNA�Ļ��򣬵��ǻ������м׻�������Ч������һ��
%�����ڽӾ�����ΪmRNA�Ļ���Ӳ��컯�Ľ������
%Dif_methylation=[W(2:end,1) num2cell(V)];
%Dif_methylation=sortrows(Dif_methylation,1);
%dmggene=Dif_methylation(:,1);
%�׻�������������ks��֤
%ɾ���׻������ݱ�ת¼�����ݶ�Ļ����
methyl_gene=met_gene;
methyl_gene=intersect(methyl_gene,gene);
methyl_gene=sortrows(methyl_gene,1);
% ɾ�����ڵĻ����
location=0;
for j=1:2 %��һ�к͵ڶ��ж���һ��
    for i=1:size(pair_co,1)
    %ת¼���еĻ����������methy1_gene�У�ɾ��
    index=find(strcmp(pair_co{i,j},methyl_gene)==1); %��i���еĻ����ڼ׻�����indexλ��
    if isempty(index)
        loc=i;
        location=[location;i]; %��¼ɾ��λ��
    end
    end
end
location(1,:)=[];
pair_co(location,:)=[];
% mRNAgene�����׻��������Ӧ
ex=ones(length(gene),1);
d=0.85;
Pair_co=pair_co;
[adj_m,weight]=create_adjacency_matrix_weight(gene,Pair_co);
adj_m(logical(eye(size(adj_m))))=0;
r=Generank(adj_m,ex,d);
r_co=num2cell(r);
relt_co=[mRNA_gene r_co];
[h_coME,p_coME]=rank_test(relt_co,gene_labs);%����
relt_co=sortrows(relt_co,-2);%���н�������
rank=(1:length(relt_co))';
rank_coME=[relt_co num2cell(rank)];
rank_coME=sortrows(rank_coME,1);%���յ���������
%%
%�����������磨cmi��
load(['./output/',cancer,'/','CMI_net.mat'])
d=0.85;
Pair_cmi=pair_cmi;
[adj_m,weight]=create_adjacency_matrix_weight(gene,Pair_cmi);
adj_m(logical(eye(size(adj_m))))=0;
ex=ones(length(gene),1);
r=Generank(adj_m,ex,d);
r_cmi=num2cell(r);
relt_cmi=[gene r_cmi];
[h_cmi,p_cmi]=rank_test(relt_cmi,gene_labs);%����
relt_cmi=sortrows(relt_cmi,-2);%���н�������
rank=(1:length(relt_cmi))';
rank_cmi=[relt_cmi num2cell(rank)];
rank_cmi=sortrows(rank_cmi,1);%���յ���������%%
%%
%ceRNA����
load(['./output/',cancer,'/','ceRNA_net.mat'])
d=0.85;
Pair_ceRNA=pair_ceRNA;
[adj_m,weight]=create_adjacency_matrix_weight(gene,Pair_ceRNA);
%�뽫��ʼֵ��Ϊ���ֵ��Ȩ��֮��
ex=ones(length(gene),1);
adj_m(logical(eye(size(adj_m))))=0;
r=Generank(adj_m,ex,d);
r_ceRNA=num2cell(r);
relt_ceRNA=[gene r_ceRNA];
[h_ceRNA,p_ceRNA]=rank_test(relt_ceRNA,gene_labs);%����
relt_ceRNA=sortrows(relt_ceRNA,-2);%���н�������
rank=(1:length(relt_ceRNA))';
rank_ceRNA=[relt_ceRNA num2cell(rank)];
rank_ceRNA=sortrows(rank_ceRNA,1);%���յ���������
%%
load(['./data/',cancer,'/',cancer,'_somatic_mutation.mat']);%miRNA����
data=zeros(length(mRNA_gene),size(sm_patient_data,2));
sm2mRNA_index=[];
for i=1:length(gene)
    index=find(strcmp(sm_gene(:,1),gene{i}));
    if index
    data(i,:)=sm_patient_data(index,:);
    end
end
A=data;
B =A';
C = A * B;
normA = sqrt(sum(A .^ 2, 2)); 
normB = sqrt(sum(B .^ 2, 1)); 
C = bsxfun(@rdivide, bsxfun(@rdivide, A * B, normA), normB);

C=C-diag(diag(C));%��Ϊ�Խ���Ϊ��ľ���

%�õ��ڽӾ���
adj=C;
adj(isnan(adj))=0;
C=zeros(length(adj),length(adj));
adj_new=zeros(length(adj),length(adj));
n=10;%ѡȡǰn�������
for i=1:length(adj)
    temp_val=adj(:,i);
    [B,index] = sortrows(temp_val,-1);
    adj_new(index(B(1:n)~=0),i)=1;
    if index(B(1:n)~=0)
    C(index(B(1:n)~=0),i)=B(1:length(index(B(1:n)~=0)));
    end
end
ex=ones(length(gene),1);
d=0.85;
r=Generank(adj_new,ex,d);%�õ�generankֵ
r_sm=num2cell(r);
relt_sm=[gene r_sm];
[h_sm,p_sm]=rank_test(relt_sm,gene_labs);%����
relt_sm=sortrows(relt_sm,-2);%���н�������
rank=(1:length(relt_sm))';
rank_sm=[relt_sm num2cell(rank)];
rank_sm=sortrows(rank_sm,1);%���յ���������

h_f=[h_coTR,h_coME,h_cmi,h_ceRNA,h_sm];
p_f=[p_coTR,p_coME,p_cmi,p_ceRNA,p_sm];
rank_f=[rank_coTR(:,3),rank_coME(:,3),rank_cmi(:,3),rank_ceRNA(:,3),rank_sm(:,3)];
rank_i=[rank_coTR(:,2),rank_coME(:,2),rank_cmi(:,2),rank_ceRNA(:,2),rank_sm(:,2)];
save(['./output/',cancer,'/',cancer,'_rank(cgc).mat'],'h_f','p_f','rank_f','gene_labs','rank_i','-v7.3'); 
%xlswrite('sy_rank.xlsx',rank_f,'Sheet1','A1');
