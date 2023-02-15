function dome_sy(cancer)
%Ԥ��ʱ��30h
addpath(genpath('mRMR_0.9_compiled2'));%��Ӽ����غ��������ļ���·��
%load GSE2034_ma2;
% load Background_vector;
% load('BLCA_gene_expression_RNAseq.mat');%ת¼������
% load('BLCA_clinic.mat');%�ٴ�����
%%
%ѡ����Ч�Ĳ�������
%data=ma2(:,co>=0);
%co=co(co>=0); 
[dmfs_e dmfs_time sample_ID ID datad, proSet]=loda_data(cancer);
%mapping�������Ա�����Ԥ�������֪������
[data ,phenotype]=mapping(dmfs_e,dmfs_time,sample_ID,ID,datad);
co=phenotype;
ma2_geneName=proSet;

[u,u1,u2,u3,u4]=background_vector_calculation(data,co);%���ɱ����ֲ�

%�ֱ�������ɸѡ�����Ĳ��˵�����
gene_num=length(data(:,1));
patient_num=length(data(1,:));
%%

tic
medians=median(data,2);%���ma2ÿһ�е���ֵ

for i=1:gene_num
    ma(i,:)=imbinarize(data(i,:),medians(i,1));%ʹ���ݶ�ֵ��
end
ma=double(ma);
%%
%n��������3���������ж��ٻ����
n1=0;
n2=0;
n3=0;
%�����п��ܵĻ���Խ��м���
for i=1:gene_num-1
    for j=i+1:gene_num 
        A=ma(i,:);
        B=ma(j,:);
        for k=1:patient_num
            if A(k)==1&B(k)==1
                C1(k)=1;  
            else
                C1(k)=0;
            end
            if A(k)==1|B(k)==1
                C2(k)=1;
            else
                C2(k)=0;
            end
            if A(k)==1&B(k)~=1
                C3(k)=1;
            else
                C3(k)=0;
            end
            if A(k)~=1&B(k)==1
                C4(k)=1;
            else
                C4(k)=0;
            end
        end 
%         C1(find((A+B)==2))=1;
%         C2(find((A+B)>=1))=1;
%         C3(find(A==1))=1;
%         C3(find(B==1))=0;
%         C4(find(A==0))=1;
%         C4(find(B==0))=0;
        U_a=(entropy(co')+entropy(A)-jointentropy(co',A))/entropy(co');
        p_a=sum(U_a<=u)/10000;
        U_b=(entropy(co')+entropy(B)-jointentropy(co',B))/entropy(co');
        p_b=sum(U_b<=u)/10000;
        U_c1=(entropy(co')+entropy(C1)-jointentropy(co',C1))/entropy(co');
        U_c2=(entropy(co')+entropy(C2)-jointentropy(co',C2))/entropy(co');
        U_c3=(entropy(co')+entropy(C3)-jointentropy(co',C3))/entropy(co');
        U_c4=(entropy(co')+entropy(C4)-jointentropy(co',C4))/entropy(co');
        p_c1=sum(U_c1<=u1)/10000;
        p_c2=sum(U_c2<=u2)/10000;
        p_c3=sum(U_c3<=u3)/10000;
        p_c4=sum(U_c4<=u4)/10000;
        %relation������������ı�
        if p_a>0.001&p_b>0.001
            if p_c1<=0.001
                n1=n1+1;
                relation1_gene{n1,1}=ma2_geneName{i,1};
                relation1_gene{n1,2}=ma2_geneName{j,1};
                relation1_gene{n1,3}=p_c1;
            end
             if p_c2<=0.001
                n2=n2+1;
                relation2_gene{n2,1}=ma2_geneName{i,1};
                relation2_gene{n2,2}=ma2_geneName{j,1};
                relation2_gene{n2,3}=p_c2;
             end
             if p_c3<=0.001|p_c4<=0.001
                n3=n3+1;
                relation3_gene{n3,1}=ma2_geneName{i,1};
                relation3_gene{n3,2}=ma2_geneName{j,1};
                relation3_gene{n3,3}=p_c3;
             end
        end
    end
end
toc
save(['./output/',cancer,'/synergistic.mat'], 'relation1_gene','relation2_gene','relation3_gene','-v7.3') 
end

% fid=fopen('syf.txt','wt');
% for i=1:size(relation3_gene,1)%��
%     fprintf(fid,'%s ','1');
%     for j=1:size(relation3_gene,2)%��
%     a = relation3_gene(i,j);
%     a = cell2mat(a);
%     fprintf(fid,'%s ',a);
%     end
%     b = relation3_val(i);
%     fprintf(fid,'%.3f',b);
%     fprintf(fid,'\n');%�ӻ��з�
% end
