function [G,valid]=H2G(H)
%  H=[1 1 1 0 0 0 0 0;0 0 0 1 1 1 0 0;1 0 0 1 0 0 1 0;0 1 0 0 1 0 0 1];
[m,n]=size(H);
valid=1;
for k=1:m                                %逐行进行高斯消元,使前m行×m列形成单位阵,从而使校验矩阵写成[I|P]形式
    vec=[k:n];                  
    if (H(k,k)==0)                       %高斯消元使H(k,k)==1
        a=find(H(k+1:m,k)~=0);
        if isempty(a)
            valid=0;
            break
        end
        a_major=a(1);
        x=k+a_major;
        H(k,vec)=rem(H(x,vec)+H(k,vec),2);
    end
    a=find(H(:,k)~=0)';                  %高斯消元使的第k列初H(k,k)==1外其余位置为0
    for x=a
          if x~=k
             H(x,vec)=rem(H(x,vec)+H(k,vec),2);
          end
    end
end
P=H(:,m+1:n);
I=eye(n-m);
G=cat(2,P',I);                           %[P'|I]即为生成矩阵
%t=mod(G'*H,2);      