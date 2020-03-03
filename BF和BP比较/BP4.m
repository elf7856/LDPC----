function vHat = BP4(rx, H,sigma, iteration)
%  rx        : Received signal vector (column vector)接收的行向量
%  H         : LDPC matrix
%  N0        : Noise variance  噪声方差
%  iteration : Number of iteration
%  vHat      : Decoded vector (0/1) 
[M N] = size(H);
% Prior probabilities先验概率
% a=ones(size(rx));
% N0=sigma^2;
   P0 = (1 + exp(-2*(sigma^(-2)).*rx)).^(-1);
   P1 = (1 + exp(2*(sigma^(-2)).*rx)).^(-1);
% P1 = ones(size(rx))./(1 + exp(-2*rx./(N0/2)));
% P0 = 1 - P1;
% Initialization初始化
K0 = zeros(M, N);
K1 = zeros(M, N);
rji0 = zeros(M, N);
rji1 = zeros(M, N);
a=repmat(P0, M, 1);
qij0 = H.*repmat(P0, M, 1);%得到初始的概率矩阵，元素为0或P0
qij1 = H.*repmat(P1, M, 1);%得到初始的概率矩阵，元素为0或P1


% Iteration迭代
for n = 1:iteration
   %fprintf('Iteration : %d\n', n); 
   % ----- Horizontal step 水平步骤-----
   for i = 1:M
      % Find non-zeros in the column
      c1 = find(H(i, :));
      
      for k = 1:length(c1)
         
         % Get column products of
         % drji\c1(l)得到Qij的值，Qij是除了校验式J以外其他校验式可信度信息已知下，信息比特为0/1的概率
         %求出每一行的每个非零元素的Qij。用于初始化。
         drji = 1;
         for l = 1:length(c1)
            if l~= k
               drji = drji*(qij0(i, c1(l)) - qij1(i, c1(l)));%%%%%%%%%%%%%%%%
            end
         end % for l
         
         rji0(i, c1(k)) = (1 + drji)/2;%由证明算出的rij0的更新方式。drji是联合概率值（乘积形式），由上面的循环得到。
         rji1(i, c1(k)) = (1 - drji)/2;
         
      end % for k
      
   end % for i以上应该是将rij初始化完毕
   
   % ------ Vertical step ------
   for j = 1:N
      
      % Find non-zeros in the row
      r1 = find(H(:, j));
      
      for k = 1:length(r1)
        
         % Get row products of prod Of
         % rij\ri(l)%rij(0/1)是假设ti=0/1下，其他参与该校验式的信息比特概率分布已知条件下，校验式满足的概率。
         %应该是除了第j位，其他元素的联合概率值（乘积形式）。
         prodOfrij0 = 1;
         prodOfrij1 = 1;   
         for l = 1:length(r1)
            if l~= k
               prodOfrij0 = prodOfrij0*rji0(r1(l), j);
               prodOfrij1 = prodOfrij1*rji1(r1(l), j);
            end
         end % for 1
         % Update constants每次更新P(a)*rij(a)，其他信息节点已知条件下，信息比特为0/1，校验式满足的概率
         K0(r1(k), j) = P0(j)*prodOfrij0;%信息比特为0,校验式满足的概率
         K1(r1(k), j) = P1(j)*prodOfrij1;%信息比特为1，校验式满足的概率
         
         % Update qij0 and qij1其他校验式已知下，信息比特为0/1的概率
         qij0(r1(k), j) = K0(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));
         qij1(r1(k), j) = K1(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));
               
      end % for k第J行的k个元素，即在每个校验式里，每个信息比特为0/1的概率
      
      % Update constants
      Ki0 = P0(j)*prod(rji0(r1, j));%每列的成率乘积，即每个信息比特参加的校验式的联合概率值
      Ki1 = P1(j)*prod(rji1(r1, j));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      % Get Qj得到最后的信息比特为0/1的概率
      Qi0 = Ki0/(Ki0 + Ki1);
      Qi1 = Ki1/(Ki0 + Ki1);
      
      % Decode Qj判决        
      if Qi1 > Qi0
         vHat(j) = 1;
      else
         vHat(j) = 0;
      end
         
   end % for j
   
end % for n

