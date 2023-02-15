%����������Ԥ��15mins
%load GSE2034_ma2;
%�׻������ݻ��߻���������
%�޸���һЩ����ĵڶ��汾
[gene,data]=loda_data_co(cancer)
ma2=data;
ma2_geneName=gene;
zee=find(ma2(:,1)==0);
azee=find(all(ma2==0,2));
n=length(ma2_geneName);
index=(1:n)';   %n�������λ��
tic
num=1;%���ڼ�¼����Ļ���
for i=1:n
    [co,p]=corr(ma2(i,:)',ma2'); %��i��������������л��������ϵ��
    res=[index co' p'];
    res=sortrows(res,-2); %�������ϵ���Ӵ�С����
    del=find(isnan(res(:,2))==1);%�ҵ�NAN��λ��
    if length(del)==length(res)%�������ȫΪ��������
        continue;
    end
    res(del,:)=[];
    %pair(((i-1)*5+1):i*5,1)=ma2_geneId(i,1);  %��ǰ�Ļ����ID
    pair_gene(((num-1)*5+1):num*5,1)={ma2_geneName{i,1}}; 
    %pair(((i-1)*5+1):i*5,2)=ma2_geneId(res(2:6,1),1); %����ص���������ID(��Ϊ����ص����Լ��������ǵ�2����������)
    pair_gene(((num-1)*5+1):num*5,2)={ma2_geneName{res(2:6,1),1}};
    pair_val(((num-1)*5+1):num*5,1:2)=res(2:6,2:3);  %%��5�������뵱ǰ��������ϵ����p-value 
    num=num+1;
end;
% xlswrite('co_net.xlsx',pair_gene,'Sheet1','A1');
% xlswrite('co_net.xlsx',pair_val,'Sheet1','C1');  %�����������洢��excel�ļ���
pair_co=[pair_gene,num2cell(pair_val)];

save(['./output/',cancer,'/co_net.mat'], 'pair_co','-v7.3') 
toc;
