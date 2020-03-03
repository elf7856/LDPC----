%利用Mackay的构造法1A构造校验矩阵
function [H]=getH(m,n)
row_flag(1:m)=0;
h=zeros(m,n);
bits_per_col=3;                                                            %使每列随机产生3个1即列重为3
for i=1:n
    a=randperm(m);
    for j=1:bits_per_col
        h(a(j),i)=1;
      row_flag(a(j))= row_flag(a(j))+1;
    end
end
max_ones_per_row=ceil(n*bits_per_col/m);                                   %求每行1的最多个数即行重的最大值
for i=1:m
    if row_flag(i)==0
        for k=1:2
            j=unidrnd(m);
            while  h(i,j)==1
                j=unidrnd(m);
            end
            h(i,j)=1;
            row_flag(i)=row_flag(i)+1;
        end
    end
    if  row_flag(i)==1
         j=unidrnd(m);
            while  h(i,j)==1
                j=unidrnd(m);
            end
            h(i,j)=1;
            row_flag(i)=row_flag(i)+1;
    end
end
          
for i=1:m                                                                  %尝试在列上分散1的位置,使行重分布尽量均匀
    j=1;
    a=randperm(n);
    while row_flag(i)>max_ones_per_row                                     %如果该行行重大于行重的最大值,则进行处理
        if h(i,a(j))==1                                                    %随机选择某一该行上为1的列来处理,将改行上的1分散到其他行上
            newrow=unidrnd(m);                                             %随机查找该列上适合放置1(行重小于最大值且该位置为0)的行
            k=0;
            while (row_flag(newrow)>=max_ones_per_row|h(newrow,a(j))==1)&k<m
                newrow=unidrnd(m);
                k=k+1;
            end
            if h(newrow,a(j))==0                                           %将待处理行的1放到找到的行上,并对两行的行标志作相应的处理
                h(newrow,a(j))=1;
                row_flag(newrow)=row_flag(newrow)+1;
                h(i,a(j))=0;
                row_flag(i)=row_flag(i)-1;
            end
        end
        j=j+1;
    end
end

for loop=1:100                                     %常是删除短环4
    success=1;
    for r=1:m
        ones_position=find(h(r,:)==1);
        ones_count=length(ones_position);
        for i=[1:r-1,r+1:m]
            common=0;
            for j=1:ones_count
                if h(i,ones_position(j))==1
                    common= common+1;
                    if  common==1
                        a=ones_position(j);         %两列的第一个共同1元素
                    end
                end
                if common==2
                    success=0;
                    common=common-1;
                    if (round(rand)==0)            %随即决定保留前面的列还是后面的列
                        b=a;                       %保留后面的列,处理前面列上的第二个共同的1元素
                        a=ones_position(j);
                    else b=ones_position(j);       %保留前面的列,处理后面列上的第二个共同的1元素
                    end
                    h(i,b)=3;                      %将该1置为3以使在以后的尝试删除中不用该值
                    newrow=unidrnd(m);
                    iteration=0;
                    while h(newrow,b)~=0 & iteration<5 %尝试5次在待交换的列中随机查找0
                        newrow=unidrnd(m);
                        iteration=iteration+1;
                    end
                    if iteration>=5                   %超过5次则随机查找非1的0或3
                        while h(newrow,b)==1
                           newrow=unidrnd(m);
                        end
                    end
                    h(newrow,b)=1;                    %将该列中找到的0或3置为1
                end
            end
        end
    end
    if success                                        %如果本次循环已不存在短环4(seccess=1),则结束循环loop
        break
    end
end
h=h==1;                                               %用0替代剩余的3,并得到构造好的校验矩阵H
H=h;