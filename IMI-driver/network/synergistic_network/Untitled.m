

n=floor(length(relation3_gene)/1000000)
% xlswrite('syn_net.xlsx',relation3_gene(1:1000000,1:2),'Sheet1','A1');
% xlswrite('co_net.xlsx',pair_val,'Sheet1','C1');  %�����������洢��excel�ļ���

% col=1;
% for i=1:n
%     xlrange=['A',num2str(col)];
% %�洢����е�λ�ã�һ�δ�һ�У��������Nocode{i}����������������Ȼ�洢�Ǿ�ת����
%     xlswrite('syn_net.xlsx',relation3_gene(1000000-1:i*1000000,1:2),'Sheet1',xlrange);
%     col=col+1000000;
%      %�洢ÿ������
% end
% xlrange=['A',num2str(col)];
% xlswrite('syn_net.xlsx',relation3_gene(col:end,1:2),'Sheet1',xlrange);
% save sy.txt -ascii relation3_val


fid=fopen('syf.txt','wt');
for i=1:size(relation3_gene,1)%��
    fprintf(fid,'%s ','1');
    for j=1:size(relation3_gene,2)%��
    a = relation3_gene(i,j);
    a = cell2mat(a);
    fprintf(fid,'%s ',a);
    end
    b = relation3_val(i);
    fprintf(fid,'%.3f',b);
    fprintf(fid,'\n');%�ӻ��з�
end
