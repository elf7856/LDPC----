function vHat = BP3(rx, H, iteration)
% Simplified log-domain sum product algorithm LDPC decoder
%◊Ó–°∫Õ“Î¬Î
%  rx        : Received signal vector (column vector)
%  H         : LDPC matrix
%  iteration : Number of iteration
%  vHat      : Decoded vector (0/1) 
[M N] = size(H);

% Prior log-likelihood (simplified). Minus sign is used for 0/1 to -1/1 mapping
Lci = rx;

% Initialization

Lrji = zeros(M, N);
%Pibetaij = zeros(M, N);

% Asscociate the L(ci) matrix with non-zero elements of H
Lqij = H.*repmat(Lci, M, 1);

for n = 1:iteration
   
   % Get the sign and magnitude of L(qij)   
   alphaij = sign(Lqij);  
   betaij = abs(Lqij);

   % ----- Horizontal step -----
   for i = 1:M
      
      % Find non-zeros in the column
      c1 = find(H(i, :));
      
      % Get the minimum of betaij
      for k = 1:length(c1)
         
         % Minimum of betaij\c1(k)
         minOfbetaij = realmax;
         for l = 1:length(c1)
            if l ~= k  
               if betaij(i, c1(l)) < minOfbetaij
                  minOfbetaij = betaij(i, c1(l));
               end
            end           
         end % for l
                 
         % Multiplication alphaij\c1(k) (use '*' since alphaij are -1/1s)
         prodOfalphaij = prod(alphaij(i, c1))*alphaij(i, c1(k)); 
         
         % Update L(rji)
         Lrji(i, c1(k)) = prodOfalphaij*minOfbetaij;
      
      end % for k
      
   end % for i

   % ------ Vertical step ------
   for j = 1:N

      % Find non-zero in the row
      r1 = find(H(:, j));
      
      for k = 1:length(r1)         
         
         % Update L(qij) by summation of L(rij)\r1(k)
         Lqij(r1(k), j) = Lci(j) + sum(Lrji(r1, j)) - Lrji(r1(k), j);
      
      end % for k
      
      % Get L(Qij)
      LQi = Lci(j) + sum(Lrji(r1, j));
      
      % Decode L(Qi)
      if LQi < 0
         vHat(j) = 1;
      else
         vHat(j) = 0;
      end
                 
   end % for j

end % for n
