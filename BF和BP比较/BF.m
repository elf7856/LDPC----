function est_code = BF(rec,H,IterNum)      
% output
%   u - decoded message
%   ite - number of iteration

% input
%   H - parity check matrix
%   re - received word
%   max_ite - maximum iteration

tic             % start timer

%%%%%%%%%% Step 1. 形成错误图样 %%%%%%%%%%
[row,col] = size(H);
hard = zeros(1,col);
y_re = rec;
iteration = 0;

% hard decision from BPSK
% y_re > 0 --> 1
% y_re <= 0 --> 0
for i = 1:col
    if y_re(i) < 0.0
        hard(i) = 1;
    else
        hard(i) = 0;
    end 
end 

y_re = hard;
syn = mod(y_re * H',2);           % 错误图样
while (sum(sum(syn)) ~= 0) && (iteration < IterNum)  %检查 if syn=0 or 达到最大迭代次数
    iteration = iteration + 1;
    %%%%%%%%%% Step 2.  计算得到每个变量节点的校验式 %%%%%%%%%%
	S=zeros(1,1);
	for i = 1:col
        v = H(:,i);
		S(i) = syn * v;
	end 
   
    %%%%%%%%%% Step 3: 存储校验式不满足个数最多的那个变量节点用于翻转 %%%%%%%%%%
    [srow,scol] = size(S);
	bflip = 1;
	flip_count = 1;

	for i = 1:scol - 1
		if S(i+1) >= S(bflip)
			bflip = i + 1;
			flip_count = flip_count + 1;
		end 
	end 

	if S(1) == S(bflip)
		bflip = 1;
	end 
    
    %%%%%%%%%% Step 4: 翻转节点 %%%%%%%%%%
    y_re(bflip)=not(y_re(bflip));
    syn = mod(y_re*H',2);  % 重新计算错误图样，即是否满足校验式
%     a = max(bflip);
%     b = y_re(1,a);
% %      f=find(a);
%     if b == 0
%         y_re(1,a) = 1;
%     else 
%         y_re(1,a) = 0;
%     end
%       syn = mod(y_re*H',2);  % 重新计算错误图样，即是否满足校验式  
end 

    est_code = y_re;
    
% 返回结果值
if (sum(sum(syn)) == 0)
    disp('BF DECODING IS SUCCESSFUL')
    est_code = y_re;
    ite = iteration;
end 

if (sum(sum(syn)) ~= 0)
    ite = iteration;
    disp('BF DECODING IS UNSUCCESSFUL')
end 

toc             % end timer