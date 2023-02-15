function demo_ceRNA(cancer)
load(['./data/',cancer,'/',cancer,'_miRNA.mat']);%miRNA����
load(['./data/',cancer,'/',cancer,'_gene_expression_RNAseq.mat'])%load����������


data=mRNA_patient_data;
gene=mRNA_gene;
ID=mRNA_patient_ID;
miRNA_data=miRNA_patient_data;
miRNA_ID=miRNA_patient_ID;

[mirna cm]=xlsread(['./data/',cancer,'/',cancer,'miRNAid.xlsx']);
cm=cm(2:end,1);
miRNA_name=cm;%ת�����ͺ�miRNA������


% load('change_name.mat');
%cm=readtable(".\change_name.csv",'ReadVariableNames',false);
% miRNA_name=cm.Var1;%ת�����ͺ�miRNA������
[new_miR,new_ge]=new_miRNA_gene(miRNA_name, miRNA_data);%�õ�����miRNA��ת¼�����ݵ�gene��miRNA��Ӧ��ϵ
 ge_ulist=unique(new_ge);
 mi_ulist=unique(new_miR);
%%
%��������mapping
[ind_ID,ind_mi]=mapping_sample(ID,miRNA_ID);
miRNA_data=miRNA_data(:,ind_mi);
data=data(:,ind_ID);
%%
%���а�֢ת¼���Ӧ�Ļ�����һ���ģ����������Դ�Ϊ.mat�ļ�
%������Ӧ��miRNA����
[g2mi]=miRNA_set(new_miR,new_ge);
%load('g2mi.mat');%�����Ӧ��miRNA���ϣ�����TCGA,miRNA���ݵõ��ļ��ϣ�

%%
pair_nums=0;
U_len=length(mi_ulist);%ȫ������,�޸ĺ�Ϊ��ȥ����һ��ΪNA��miRNA��Ŀ

%%
% cor_gene_name={};
% 
% for i=1:length(miRNA_name)
% [co_m_g,p_m_g]=corr(miRNA_data(i,:)',data'); %��i��miRNA�ͻ��������ϵ��
% %co��col��ʾ��ÿ����������ϵ��
% %�ҳ������Ϊ��ֵ����pֵС��0.05�Ļ����
% %cor_gene_name���򼯺�
% pdco=(co_m_g<0)+(p_m_g<0.05);
% cor_gene_index=find(pdco>=2);
% cor_gene_name=[cor_gene_name;gene(find(pdco>=2))];%������
% cor_gene_name=unique(cor_gene_name);
% end
[mi2g]=miRNA_set(new_ge,new_miR);
[cor_gene_name]=miRNA_gene_corr(gene,miRNA_name,mi_ulist,mi2g,data,miRNA_data);


%%
%
candi_gene_num=length(cor_gene_name);
tic
for j=1:candi_gene_num-1
    %%
    indx_A=find(strcmp(cor_gene_name{j},ge_ulist)==1);
    if length(indx_A) %��A���ϲ�Ϊ�������B�����Ӧ��miRNA����
    set_A=[];  
    set_A=g2mi{indx_A,1};%�õ��û����Ӧ��miRNA�ļ���
        for k=j+1:candi_gene_num
            %%
            indx_B=find(strcmp(cor_gene_name{k},ge_ulist)==1);
            if length(indx_B)%���B�����Ӧ��miRNA���ϲ�Ϊ�գ�����AB�Ľ���
                set_B=[];
                set_B=g2mi{indx_B,1};%�õ��û����Ӧ��miRNA�ļ���
                inter_AB=intersect(set_A,set_B);%��AB�Ľ���
                if length(inter_AB) %��������Ϊ��,����г����ηֲ�����
                    %p=1-hygecdf(��������-1,ȫ������,����A����,����B����)��
                    p_hyp=1-hygecdf(length(inter_AB)-1,U_len,length(set_A),length(set_B));
                    if p_hyp<0.05
                        [co_g_g,p_g_g]=corr(data(j,:)',data(k,:)');
                        if co_g_g>0&&p_g_g<0.05
                            %��¼�����
                            pair_nums=pair_nums+1;
                            pair_ceRNA{pair_nums,1}=cor_gene_name{j};
                            pair_ceRNA{pair_nums,2}=cor_gene_name{k};
                            pair_ceRNA{pair_nums,3}=co_g_g;
                        end
                    end
                end
            end
        end
    end
end
 toc
save(['./output/',cancer,'/ceRNA_net.mat'],'pair_ceRNA','-v7.3'); 
end
