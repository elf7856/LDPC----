m=9;n=12;
H=getH(m,n);
[G,valid]=H2G(H); 
while valid==0                   %valid作为校验矩阵是否为满秩的标志,若为非满秩(valid=0),则返回重新利用1A随机构造
H=getH(m,n);                     %构造校验矩阵H
[G,valid]=H2G(H);                %将校验矩阵H转化为生成矩阵G
end

% t=mod(G*H,2);
% [i.j]=find(t);

save G