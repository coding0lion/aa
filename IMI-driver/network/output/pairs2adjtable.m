%%�ɹ�ϵ��ת��Ϊ�ڽӱ�(���Լ����Լ�)
function pairs2adjtable(cancer,netname,pair)
    load(['./data/',cancer,'/',cancer,'_gene_expression_RNAseq.mat'])%load����������
    gene=mRNA_gene;
    temp_a=pair(:,1:2);
    num=1;
    a={};
    for i=1:length(temp_a)
    if ~isempty(find(strcmp(temp_a{i,1}, gene)==1))&&~isempty(find(strcmp(temp_a{i,2}, gene)==1))
        a(num,:)=temp_a(i,:);
        num=num+1;
        end
    end
    % a=relation3_gene(:,1:2);
    %��Ϊ�Գƾ���
    aa=[a(:,2),a(:,1)];
    haved=unique(a);
    isolate=setdiff(gene,haved);%�õ�������
    isolate_pair=[gene,gene];
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
    e = e-eye(20530);
    adj_m=e;

    [ii, jj] = find(adj_m); % row and col indices of connections 
    y = accumarray(ii, jj-1 , [], @(x){sort(x.')}); % get all nodes connected to each node, 
    node=[0:1:length(gene)-1]';

    fid=fopen(['output/',cancer,'/','sub_',cancer,'_',netname,'.txt'],'wt');
    for i=1:size(gene,1)%��
    b = node(i);
    fprintf(fid,'%.0f ',b);
    if i<=size(y,1)
    for j=1:size(y{i},2)%��
    a = y{i}(j);
    %     a = cell2mat(a);
    fprintf(fid,'%.0f ',a);
    end
    end
    fprintf(fid,'\n');%�ӻ��з�
    end
end