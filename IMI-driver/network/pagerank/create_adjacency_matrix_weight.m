%�޸ĺ󣬿��Է���0/1�ڽӾ������Ȩ�ص��ڽӾ���
function [adj_m,adj_m_weight]=create_adjacency_matrix(gene,Pair)
% n=length(gene);
% temp=zeros(n,n);
% tic
% for i=1:length(Pair)
% col=find(strcmp(Pair{i,1},gene)==1);
% row=find(strcmp(Pair{i,2},gene)==1);
% temp(col,row)=1;
% temp(row,col)=1;
% end
% toc
% pa=Pair(1:end,1:2);
% ge=unique(pa);
% set_set=setdiff(gene,ge);
% for j=1:length(set_set)
%     x=find(strcmp(set_set{j},gene)==1);
%     temp(x,x)=1;
% end
tic
a=Pair(:,1:2);
weight=Pair(:,3);
weight=cell2mat(weight);
%��Ϊ�Գƾ���
aa=[a(:,2),a(:,1)];
haved=unique(a);
isolate=setdiff(gene,haved);%�õ�������
isolate_pair=[isolate,isolate];
final_pair=[a;aa;isolate_pair];

[~,~,k]=unique(final_pair);
a=reshape(k,[],2);
b=unique(a);
c=find(b);
%�������ֵĹ�ϵ��
for i=1:length(b)
idx=find(a==b(i));
a(idx)=c(i);
end
e=accumarray(a,1);
%%
%���ܻ���ֶ�Ӧ��ȫ������ת�ú�ά�Ȳ�ͬ�������˷����ٶȺܿ죩
% e1=e-diag(diag(e));%��Ϊ�Խ���Ϊ��ľ���
% e=e+e1';%��ȫ�Գƾ���
%%
adj_m=e;
%Ȩ�ؾ���
pair_weight=a(1:length(weight),:);
e_weight=zeros(length(gene),length(gene));
nn=(pair_weight(:,2)-1).*20530+pair_weight(:,1);
e_weight(nn)=weight;
e_weight=e_weight+e_weight';
adj_m_weight=e_weight;
toc
end