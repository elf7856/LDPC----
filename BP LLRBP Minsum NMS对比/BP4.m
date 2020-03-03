function vHat = BP4(rx, H,sigma, IterNum)
%%%概率域译码
%计算译码参数
%先根据接收数据和信道信息初始化每个比特节点的概率软信息，再将比特节点的概率信息传递给概率信息传递给校验节点，校验节点的概率信息是参与该校验方程的
%所有比特节点对该校验方程的概率贡献，如何校验节点将概率信息传递给比特节点，对于某个比特节点而言，该软信息就是他所参与的所有校验方程对他的概率贡献，
%最后，比特节点利用从他参与的所有校验方程那里得到概率信息以及信道初始化细腻译码判决，如此完成了一次信息传递和迭代过程，当某次判决译码成功了或者达
%到迭代次数上限时就停止译码。否则进入下一次迭代过程
%  rx        : 接收的含噪声数据向量
%  H         : 校验矩阵
%  sigma     : 信道噪声标准差
%  IterNum   : 迭代次数
%  vHat      : 译码结果 
   [M N] = size(H);
% 先验概率
   P0 = (1 + exp(-2*(sigma^(-2)).*rx)).^(-1);
   P1 = (1 + exp(2*(sigma^(-2)).*rx)).^(-1);
% 初始化
   K0 = zeros(M, N);
   K1 = zeros(M, N);
   rji0 = zeros(M, N);
   rji1 = zeros(M, N);
   a=repmat(P0, M, 1);
   qij0 = H.*repmat(P0, M, 1);%得到初始的概率矩阵，元素为0或P0
   qij1 = H.*repmat(P1, M, 1);%得到初始的概率矩阵，元素为0或P1
   success=0;z=0;
%迭代
while((success==0)&(z<IterNum))    
for n = 1:IterNum
   %%%% 水平步骤 %%%%
   for i = 1:M
      c1 = find(H(i, :));  % 在行找1
        for k = 1:length(c1)    
           % drji\c1(l)得到Qij的值，Qij是除了校验式J以外其他校验式可信度信息已知下，信息比特为0/1的概率
           %求出每一行的每个非零元素的Qij。用于初始化。
           drji = 1;
          for l = 1:length(c1)
            if l~= k
               drji = drji*(qij0(i, c1(l)) - qij1(i, c1(l)));
            end
          end % for l
           rji0(i, c1(k)) = (1 + drji)/2;%由证明算出的rij0的更新方式。drji是联合概率值（乘积形式），由上面的循环得到。
           rji1(i, c1(k)) = (1 - drji)/2;  
        end % for k
   end % for i以上是将rij初始化完毕
   % ------ 垂直步骤 ------
   for j = 1:N
      % 在列找1
      r1 = find(H(:, j));   
      for k = 1:length(r1)   
         % rij\ri(l)%rij(0/1)是假设ti=0/1下，其他参与该校验式的信息比特概率分布已知条件下，校验式满足的概率。
         %除了第j位，其他元素的联合概率值（乘积形式）。
         prodOfrij0 = 1;
         prodOfrij1 = 1;   
         for l = 1:length(r1)
            if l~= k
               prodOfrij0 = prodOfrij0*rji0(r1(l), j);
               prodOfrij1 = prodOfrij1*rji1(r1(l), j);
            end
         end % for 1
         % 每次更新P(a)*rij(a)，其他信息节点已知条件下，信息比特为0/1，校验式满足的概率
         K0(r1(k), j) = P0(j)*prodOfrij0;%信息比特为0,校验式满足的概率
         K1(r1(k), j) = P1(j)*prodOfrij1;%信息比特为1，校验式满足的概率
         % 更新qij0 和 qij1其他校验式已知下，信息比特为0/1的概率
         qij0(r1(k), j) = K0(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));
         qij1(r1(k), j) = K1(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));      
      end % for k第J行的k个元素，即在每个校验式里，每个信息比特为0/1的概率
      % 更新
      Ki0 = P0(j)*prod(rji0(r1, j));%每列的成率乘积，即每个信息比特参加的校验式的联合概率值
      Ki1 = P1(j)*prod(rji1(r1, j));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      %  Qj得到最后的信息比特为0/1的概率
      Qi0 = Ki0/(Ki0 + Ki1);
      Qi1 = Ki1/(Ki0 + Ki1);
      
      % 译码 Qj判决        
      if Qi1 > Qi0
         vHat(j) = 1;
      else
         vHat(j) = 0;
      end    
   end % for j 
end % for n
     if mod(vHat*H',2)==0                                  %判决译码是否成功,若成功(success=1),则退出循环运算,否则继续.
       success=1;
     break
     else  z=z+1;   
     end
end
end

