%mapping���������ݺ��ٴ�����
%�������˵��
%dmfs_e�ٴ����ݵ�Ԥ�����
%dmfs_time�ٴ����ݴ��ʱ��
%sample_ID�ٴ������е�����ID
%ID���������ݵ�����id
%ma���������ݵĸ�����ı�����������data����ʦ������ͳһ����Ϊ��ma����ʹ�ã�
function [ma ,phenotype]=mapping_met(dmfs_e,dmfs_time,sample_ID,ID,ma)
%%
%���н������������ٴ����ݵ���������ƥ��
ID=reshape(ID,length(ID),1);%�����������е�������ID
%%
%����������ݵ�Ԥ����������������1825Ϊ����
for k=1:length(sample_ID)
    if dmfs_e(k,1)==0&&dmfs_time(k,1)<1825%Ԥ��δ֪�����
        yu(k,1)=2;
    elseif dmfs_e(k,1)==1&&dmfs_time(k,1)<1825
        yu(k,1)=1;
    elseif dmfs_e(k,1)==0&&dmfs_time(k,1)>=1825
        yu(k,1)=0;
    else
        yu(k,1)=0; 
    end
end
unkonw=find(yu==2);%��¼Ԥ��δ֪������������
bad=find(yu==1);
good=find(yu==0);
%%
%ȥ�����ϸ�����ݣ�����
sample_ID(unkonw)={0};
cont=zeros(length(ID),1);%�����������������������һ��
phenotype(:,1) = dmfs_time;
phenotype(:,2) = cont;
for i=1:length(ID)
    for j=1:length(sample_ID)%�ٴ�������������ID����
        if  strcmp(ID{i,1},sample_ID{1,j})
            cont(i,1)=j;
            if find(j==bad)
                phenotype(i,2)=1;%���Ԥ�����
            end
        end
    end
end
del_sample=find(cont==0);
cont(del_sample)=[];
phenotype(del_sample,:)=[];%ɾ�����ϸ����ݺ��Ԥ�����
%%
ma(:,del_sample)=[];
end