
function [H,c]=ldpc_matrix1(msg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z=24.28.32.36.40.44.48.52.56.60.64.68.72.76.80.84.88.92.96.
% Z=24;%Rate=1/2
% rotmatrix = ...
% [-1 94 73 -1 -1 -1 -1 -1 55 83 -1 -1  7  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
% -1 27 -1 -1 -1 22 79  9 -1 -1 -1 12 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
% -1 -1 -1 24 22 81 -1 33 -1 -1 -1  0 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1;
% 61 -1 47 -1 -1 -1 -1 -1 65 25 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1;
% -1 -1 39 -1 -1 -1 84 -1 -1 41 72 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1;
% -1 -1 -1 -1 46 40 -1 82 -1 -1 -1 79  0 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1;
% -1 -1 95 53 -1 -1 -1 -1 -1 14 18 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1;
% -1 11 73 -1 -1 -1  2 -1 -1 47 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1;
% 12 -1 -1 -1 83 24 -1 43 -1 -1 -1 51 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1;
% -1 -1 -1 -1 -1 94 -1 59 -1 -1 70 72 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1;
% -1 -1  7 65 -1 -1 -1 -1 39 49 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0;
% 43 -1 -1 -1 -1 66 -1 41 -1 -1 -1 26  7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0];

% % % 
% % Z = 24; % code rate = 2/3A
% % rotmatrix = ...
% %  [3  0 -1 -1  2  0 -1  3  7 -1  1  1 -1 -1 -1 -1  1  0 -1 -1 -1 -1 -1 -1;
% % -1 -1  1 -1 36 -1 -1 34 10 -1 -1 18  2 -1  3  0 -1  0  0 -1 -1 -1 -1 -1;
% % -1 -1 12  2 -1 15 -1 40 -1  3 -1 15 -1  2 13 -1 -1 -1  0  0 -1 -1 -1 -1;
% % -1 -1 19 24 -1  3  0 -1  6 -1 17 -1 -1 -1  8 39 -1 -1 -1  0  0 -1 -1 -1;
% % 20 -1  6 -1 -1 10 29 -1 -1 28 -1 14 -1 38 -1 -1  0 -1 -1 -1  0  0 -1 -1;
% % -1 -1 10 -1 28 20 -1 -1  8 -1 36 -1  9 -1 21 45 -1 -1 -1 -1 -1  0  0 -1;
% % 35 25 -1 37 -1 21 -1 -1  5 -1 -1  0 -1  4 20 -1 -1 -1 -1 -1 -1 -1  0  0;
% % -1  6  6 -1 -1 -1  4 -1 14 30 -1  3 36 -1 14 -1  1 -1 -1 -1 -1 -1 -1  0];

Z = 24; % code rate = 2/3B
rotmatrix = ...
[2 -1 19 -1 47 -1 48 -1 36 -1 82 -1 47 -1 15 -1 95  0 -1 -1 -1 -1 -1 -1;
-1 69 -1 88 -1 33 -1  3 -1 16 -1 37 -1 40 -1 48 -1  0  0 -1 -1 -1 -1 -1;
10 -1 86 -1 62 -1 28 -1 85 -1 16 -1 34 -1 73 -1 -1 -1  0  0 -1 -1 -1 -1;
-1 28 -1 32 -1 81 -1 27 -1 88 -1  5 -1 56 -1 37 -1 -1 -1  0  0 -1 -1 -1;
23 -1 29 -1 15 -1 30 -1 66 -1 24 -1 50 -1 62 -1 -1 -1 -1 -1  0  0 -1 -1;
-1 30 -1 65 -1 54 -1 14 -1  0 -1 30 -1 74 -1  0 -1 -1 -1 -1 -1  0  0 -1;
32 -1  0 -1 15 -1 56 -1 85 -1  5 -1  6 -1 52 -1  0 -1 -1 -1 -1 -1  0  0;
-1  0 -1 47 -1 13 -1 61 -1 84 -1 55 -1 78 -1 41 95 -1 -1 -1 -1 -1 -1  0];
 
% Z = 24; 
% % % code rate = 3/4A
% rotmatrix = ...
% [6  38  3 93 -1 -1 -1 30 70 -1 86 -1 37 38  4 11 -1 46 48  0 -1 -1 -1 -1; 
% 62 94 19 84 -1 92 78 -1 15 -1 -1 92 -1 45 24 32 30 -1 -1  0  0 -1 -1 -1; 
% 71 -1 55 -1 12 66 45 79 -1 78 -1 -1 10 -1 22 55 70 82 -1 -1  0  0 -1 -1; 
% 38 61 -1 66  9 73 47 64 -1 39 61 43 -1 -1 -1 -1 95 32  0 -1 -1  0  0 -1; 
% -1 -1 -1 -1 32 52 55 80 95 22  6 51 24 90 44 20 -1 -1 -1 -1 -1 -1  0  0; 
% -1 63 31 88 20 -1 -1 -1  6 40 56 16 71 53 -1 -1 27 26 48 -1 -1 -1 -1  0] 
 
% Z = 30; % code rate = 3/4B
% rotmatrix = ...
% [-1 81 -1 28 -1 -1 14 25 17 -1 -1 85 29 52 78 95 22 92  0  0 -1 -1 -1 -1;
% 42 -1 14 68 32 -1 -1 -1 -1 70 43 11 36 40 33 57 38 24 -1  0  0 -1 -1 -1;
% -1 -1 20 -1 -1 63 39 -1 70 67 -1 38  4 72 47 29 60 5  80 -1  0  0 -1 -1;
% 64  2 -1 -1 63 -1 -1  3 51 -1 81 15 94  9 85 36 14 19 -1 -1 -1  0  0 -1;
% -1 53 60 80 -1 26 75 -1 -1 -1 -1 86 77  1  3 72 60 25 -1 -1 -1 -1  0  0;
% 77 -1 -1 -1 15 28 -1 35 -1 72 30 68 85 84 26 64 11 89  0 -1 -1 -1 -1  0]

%Z=24;% rate = 5/6
% rotmatrix = ...
% [1 25 55 -1 47  4 -1 91 84  8 86 52 82 33  5  0 36  20  4 77 80  0 -1 -1;
% -1 6 -1  36 40 47 12 79 47 -1 41 21 12 71 14 72  0  44 49  0  0  0  0 -1;
% 51 81 83  4 67 -1 21 -1 31 24 91 61 81  9 86 78 60  88 67 15 -1 -1  0  0; 
% 68 -1 50 15 -1 36 13 10 11 20 53 90 29 92 57 30 84  92 11 66 80 -1 -1  0]
H = zeros(size(rotmatrix)*Z);
Zh = diag(ones(1,Z),0);

for r=1:size(rotmatrix,1)
    for c=1:size(rotmatrix,2)
        rotidx = rotmatrix(r,c);
        if (rotidx > -1)            
            Zt = circshift(Zh,[0 rotidx]);
        else
            Zt = zeros(Z);
        end
        limR = (r-1)*Z+1:r*Z;
        limC = (c-1)*Z+1:c*Z;
        H(limR,limC) = Zt;
    end
end

[n,m]=size(H);
%计算列重
col_flag(1:m)=0;
for j=1:m
  ind=find(H(:,j)==1);
  col_flag(j)=length(ind);
end
%计算行重
row_flag(1:n)=0;
for i=1:n
  ind=find(H(i,:)==1);
  row_flag(i)=length(ind);
end
mLen = size(H,1);
cLen = size(H,2);
% Rename variables to match those used by the author
n = cLen;
%m = cLen - mLen;
m = mLen;
s = msg;
% Let c = [s p1 p2] be the codeword

Hrow1 = H(:,end);

for i=1:cLen
    if Hrow1(i) == 1
        g = i;
        break;
    end
end
g = mLen - g;

wA = n-m;
wB = g;
eA = wA;
eB = wA + wB;

% Extract the submatrices A, B, C, D, E and T
A = H(1:m-g,1:eA) ;
B = H(1:m-g,eA+1:eB);
T = H(1:m-g,eB+1:end);
C = H(m-g+1:end,1:eA);
D = H(m-g+1:end,eA+1:eB);
E = H(m-g+1:end,eB+1:end);

% Calculate p1 and p2
invT = (inv(T)); 
ET1 = -(E*invT);

Iup = diag(ones(1,size(ET1,2)),0);
Idn = diag(ones(1,size(ET1,1)),0);
X = [Iup zeros(size(Iup,1),size(Idn,2)); ET1 Idn];
Y = X*H;
spy(Y); 

phi = ET1*B + D;
xtra = ET1*A + C;
p1 = mod(phi*xtra*(s'),2)';
p2 = mod(invT*(A*(s') + B*(p1')),2)';
c = [s p1 p2];

% Checking c*H'=0;
zero = mod(c*H',2);
if sum(zero)== 0
     disp('OK!')
else
    disp('Error')
end
end

