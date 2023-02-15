
load('BLCA_gene_expression_RNAseq.mat');

gene=sort(gene);

[data,str]=xlsread('pp.xlsx');
pathway=zeros(length(gene),size(str,1));
for i=1:length(gene)
    for j=1:size(str,1)
    if find(strcmp(gene{i},str(j,:))==1)
        pathway(i,j)=1;
    end
    end
end
%����
AA=pathway(85,:);

BB =pathway(479,:);

dot(AA,BB)/(sqrt(sum(AA.*AA))*sqrt(sum(BB.*BB)));
%ͨ����������������ƶ�
%A��ά��M��N����B��N��P����ʵ���ϣ������������ļ��� - A�е���������B�е���������
%���ÿ��a��b���������ƶȵ÷֣�����a�����Ծ���A���������У�����b�������У��Ӿ���B��ʼ��
A=pathway;
B =A';
C = A * B;
normA = sqrt(sum(A .^ 2, 2)); 
normB = sqrt(sum(B .^ 2, 1)); 
C = bsxfun(@rdivide, bsxfun(@rdivide, A * B, normA), normB);
C=C-diag(diag(C));%��Ϊ�Խ���Ϊ��ľ���

%�õ��ڽӾ���
adj=C;
adj(isnan(adj))=0;
adj(adj>0)=1;

[col,row]=find(C>0);
%��ϵ��
val=C(find(C>0));
pairs=[row,col,val];
unique(pairs)

[ii, jj] = find(adj); % row and col indices of connections 
y = accumarray(ii, jj-1 , [], @(x){sort(x.')}); % get all nodes connected to each node, 
node=[0:1:length(gene)-1]';

fid=fopen('sub_pathway.txt','wt');
for i=1:size(gene,1)%��
    b = node(i);
    fprintf(fid,'%.0f ',b);
    for j=1:size(y{i},2)%��
    a = y{i}(j);
%     a = cell2mat(a);
    fprintf(fid,'%.0f ',a);
    end
    fprintf(fid,'\n');%�ӻ��з�
end