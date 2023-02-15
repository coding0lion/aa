function pair = GR_CMI(ma,gene,phenotype,permute_type,p_cut,candidate_pair)
%%%%%%%%%%�ó������ڼ������A����͵�����ԶԻ���B�������ԣ����õķ�������������Ϣ
%%%%%%��������������£�
%%%%%%ma:m*n:ת¼��������ݣ�һ��Ϊһ������һ��Ϊһ������
%%%%%%gene:m*1��ÿһ�л����Ӧ��gene symbol
%%%%%%phenotype��n*1��n�������ı���(0����1)
%%%%%permute_type:����������ϵ�Ե�ʱ��pֵ�ļ��㷽����1 (ֱ�Ӳ����õ������ֲ���Ȼ����Ƴ�ÿ�Ի���Ļ�������ֵ��pֵ) or 2 (��ÿ�Ի���ʹ��permutation�ķ��������㵱ǰ��һ�Ի���Ե�pֵ).
%%%%%��permute_type=1ʱ�������ȽϿ죨������ô׼ȷ������permute_type=2ʱ�������Ƚ�������׼ȷ��
%%%%%p_cut:p-value����ֵ�����ھ�����ǰ�Ļ����ϵ���Ƿ�����
%%%%candidate_pair�������ĺ�ѡ�����ϵ�ԣ����ò���ΪĬ��ֵʱ�����Ըò������и�ֵ������������ǰ���m����������п��ܵĹ�ϵ��
%%%����candidate_pair�����˽�����PPI��ϵ����Ϊ���룬�������Լ��ټ�������ͬʱ����PPI��ϵ�ԣ���һ���̶ȱ�������Խ��
%%%��ֻ����PPI���н����Ĺ�ϵ�ԣ�A-B���ż���A��B�Լ�B��A��������ϵ

%%�����Pair: num*4��num�������Ļ���������ϵ����һ�к͵ڶ���Ϊ����ԣ�������ΪCMIֵ��������Ϊpֵ
%%�������ǵڶ��л���Ա��͵Ļ���Ϣ�������Ե�һ�еĻ���

%%ͬʱ����Ҳ�����һ��CMI_net.mat�ļ������ļ�Ҳ�洢��Pair������

if (nargin<5)
    error('There should be at least five input parameters');
end

n=length(gene);

if (nargin==5)
    cou=0;
    for i=1:n
        for j=(i+1):n
            cou=cou+1;
            candidate_pair{cou,1}=gene{i,1};
            candidate_pair{cou,2}=gene{j,1};
            loc(cou,1)=i;
            loc(cou,2)=j;
        end
    end
else
    num=length(candidate_pair(:,1));
    flag=ones(num,1);
    loc=zeros(num,2);
    index=(1:length(gene))';
    for i=1:num
        lo=strcmpi(gene,candidate_pair{i,1});
        if (sum(lo)==1)
            loc(i,1)=index(lo,1);
            lo=strcmpi(gene,candidate_pair{i,2});
            if (sum(lo)==1)
                loc(i,2)=index(lo,1);
            else
                flag(i,1)=0;
            end
        else
            flag(i,1)=0;
        end
    end
    candidate_pair=candidate_pair(flag==1,:);
    loc=loc(flag==1,:);
end


path(path,'.\MI_tool');

if (permute_type==1)
    w=length(phenotype);
    s=fix(0.35*w);
    e=w-s+1;
    n=length(gene);
    
    OV_CMI=zeros(100000,1);
    for i=1:100000
        index=randperm(n)';
        xx=(ma(index(1,1),:))';
        yy=(ma(index(2,1),:))';
        
        yy=yy-median(yy);
        yy(yy<=0)=0;
        yy(yy>0)=1;
        
        tri=[xx yy phenotype];
        tri=sortrows(tri,1);
        
        
        h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
        h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
        
        h=h2-h1;  %CMI
        OV_CMI(i,1)=h;
    end
    pd = fitdist(OV_CMI,'Normal');
    
    n=length(candidate_pair);
    cou=0;
    for i=1:n
        xx=(ma(loc(i,1),:))';
        yy=(ma(loc(i,2),:))';
        yy=yy-median(yy);
        yy(yy<=0)=0;
        yy(yy>0)=1;
        
        tri=[xx yy phenotype];
        tri=sortrows(tri,1);
        
        h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
        h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
        
        h=h2-h1;
        
        p=(1-cdf(pd,abs(h)))*2;
        if (p<=p_cut)
            cou=cou+1;
            pair{cou,1}=candidate_pair{i,1};
            pair{cou,2}=candidate_pair{i,2};
            pair{cou,3}=h;
            pair{cou,4}=p;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        xx=(ma(loc(i,2),:))';
        yy=(ma(loc(i,1),:))';
        yy=yy-median(yy);
        yy(yy<=0)=0;
        yy(yy>0)=1;
        
        tri=[xx yy phenotype];
        tri=sortrows(tri,1);
        
        h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
        h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
        
        h=h2-h1;
        
        p=(1-cdf(pd,abs(h)))*2;
        if (p<=p_cut)
            cou=cou+1;
            pair{cou,1}=candidate_pair{i,2};
            pair{cou,2}=candidate_pair{i,1};
            pair{cou,3}=h;
            pair{cou,4}=p;
        end
    end
else
    w=length(phenotype);
    s=fix(0.35*w);
    e=w-s+1;
    n=length(candidate_pair);
    cou=0;
    for i=1:n
        xx=(ma(loc(i,1),:))';
        yy=(ma(loc(i,2),:))';
        yy=yy-median(yy);
        yy(yy<=0)=0;
        yy(yy>0)=1;
        
        tri=[xx yy phenotype];
        tri=sortrows(tri,1);
        
        h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
        h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
        
        h=h2-h1;
        
        null_cou=0;
        for j=1:1000
            lo=randperm(w)';
            tri=[xx(lo,1) yy phenotype];
            tri=sortrows(tri,1);
            
            h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
            h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
            
            null_h=h1-h2;
            if (abs(null_h)>=abs(h))
                null_cou=null_cou+1;
            end
        end
        
        p=null_cou/1000;
        
        if (p<=p_cut)
            cou=cou+1;
            pair{cou,1}=candidate_pair{i,1};
            pair{cou,2}=candidate_pair{i,2};
            pair{cou,3}=h;
            pair{cou,4}=p;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        xx=(ma(loc(i,2),:))';
        yy=(ma(loc(i,1),:))';
        yy=yy-median(yy);
        yy(yy<=0)=0;
        yy(yy>0)=1;
        
        tri=[xx yy phenotype];
        tri=sortrows(tri,1);
        
        h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
        h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
        
        h=h2-h1;
        
        null_cou=0;
        for j=1:1000
            lo=randperm(w)';
            tri=[xx(lo,1) yy phenotype];
            tri=sortrows(tri,1);
            
            h1=mutualInformation(tri(1:s,2),tri(1:s,3)); %?��35%????????
            h2=mutualInformation(tri(e:w,2),tri(e:w,3)); %?��35%????????
            
            null_h=h1-h2;
            if (abs(null_h)>=abs(h))
                null_cou=null_cou+1;
            end
        end
        
        p=null_cou/1000;
        
        if (p<=p_cut)
            cou=cou+1;
            pair{cou,1}=candidate_pair{i,2};
            pair{cou,2}=candidate_pair{i,1};
            pair{cou,3}=h;
            pair{cou,4}=p;
        end
    end
end

end



