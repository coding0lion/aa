function [cor_gene_name]=miRNA_gene_corr(gene,miRNA_name,mi_ulist,mi2g,data,miRNA_data)
%miRNA_name�Ѿ�ɾ����NA����һ���
tic
cor_gene_name={};
n=length(mi_ulist);
    for k=1:n
        index_mi=find(strcmp(mi_ulist{k},miRNA_name)==1);
        mi_data=miRNA_data(index_mi,:)';%��ȡ��miRNA����
        del_index=find(mi_data==-1);%ɾ��NA������
        mi_data(del_index,:)=[];
        if ~isempty(mi_data)
        if ~isempty(mi2g{k,1})
            index_ge=[];
            for j=1:length(mi2g{k,1})
                index_ge=[index_ge;find(strcmp(mi2g{k,1}(j),gene)==1)];%miRNA��Ӧ�Ļ�����ת¼�������ϵ�����
            end
            gene_data=data(index_ge,:)';%��ȡmiRNA��Ӧ�����ת¼������
            gene_data(del_index,:)=[];%ɾ��miRNAΪNA��Ӧλ�õ�����
            [co_m_g,p_m_g]=corr(mi_data,gene_data,'Type','Spearman');%���������
            pdco=(co_m_g<0)+(p_m_g<0.05);%�������pС��0.05
            cor_gene_index=find(pdco>=2);
            cor_gene_name=[cor_gene_name;gene(find(pdco>=2))];%������
            cor_gene_name=unique(cor_gene_name);
        end
        end
    end
toc     
end  