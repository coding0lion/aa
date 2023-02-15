

n=floor(length(relation3_gene)/1000000)
% xlswrite('syn_net.xlsx',relation3_gene(1:1000000,1:2),'Sheet1','A1');
% xlswrite('co_net.xlsx',pair_val,'Sheet1','C1');  %将共表达网络存储到excel文件中

% col=1;
% for i=1:n
%     xlrange=['A',num2str(col)];
% %存储表格中的位置，一次存一行，所以你的Nocode{i}必须是行向量，不然存储是就转下置
%     xlswrite('syn_net.xlsx',relation3_gene(1000000-1:i*1000000,1:2),'Sheet1',xlrange);
%     col=col+1000000;
%      %存储每个数据
% end
% xlrange=['A',num2str(col)];
% xlswrite('syn_net.xlsx',relation3_gene(col:end,1:2),'Sheet1',xlrange);
% save sy.txt -ascii relation3_val


fid=fopen('syf.txt','wt');
for i=1:size(relation3_gene,1)%行
    fprintf(fid,'%s ','1');
    for j=1:size(relation3_gene,2)%列
    a = relation3_gene(i,j);
    a = cell2mat(a);
    fprintf(fid,'%s ',a);
    end
    b = relation3_val(i);
    fprintf(fid,'%.3f',b);
    fprintf(fid,'\n');%加换行符
end
