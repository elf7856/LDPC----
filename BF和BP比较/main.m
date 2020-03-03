close all
clear all
%基本参数
n=576;                                                                   %%%总码长
k=288;                                                                    %%%信息位长度
% n = 1248;
% k = 576;
rate=k/n;                                                                 %%%码率
IterNum=5;

msg = round(rand(1,k));
save msg;
% [H,c]=ldpc_matrix(msg);
[H,c]=bianma(msg);
ferrlim=5;
Npf = ferrlim*n;
Ndb=12;
EbN0db=0:1:12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN1=1:length(EbN0db)
    en1=10^(EbN0db(nEN1)/10);
    sigma=1/sqrt(2*rate*en1);
    nframe1=0;
    Err1=0;
    while nframe1 < ferrlim
        nframe1=nframe1+1;
        code=c;
        I=1-2*code;
        rec1=I+sigma*randn(1,n);        
%        [est_code1,success]= LDPC_decoder(rec1,sigma,H,IterNum);  
        est_code22=BP4(rec1,H,sigma, IterNum);%%%%%%%概率域BP
        est_code2= est_code22(1:k);
        err1=length(find(est_code2~=msg));
        Err1=Err1+err1;
    end
    BER1(nEN1)=Err1/(ferrlim*k);
    if BER1(nEN1) < 1/(ferrlim*k);
        BER1(nEN1) = .1/(ferrlim*k);
    end  
end

  for nEN2=1:length(EbN0db)
    en2=10^(EbN0db(nEN2)/10);
    sigma=1/sqrt(2*rate*en2);
    nframe1=0;
    Err2=0;
    while nframe1<ferrlim
       nframe1=nframe1+1;
        code=c;
        I=1-2*code;
       rec1=I+sigma*randn(1,n);
%         I=code;
%         b=c;
%         rec=I+sigma*randn(1,n);
%         for i=1:n
%             if rec(i) < 0.5
%          rec(i) = 0;
%       else
%          rec(i) = 1;
%             end
%         end
%  b(1,5)=1;b(1,41)=1;
%  b(1,264)=1;b(1,467)=0;
%  b(1,513)=0;
%   rec=b;%+sigma*randn(1,n); 
% est_code = BP3(b, H,IterNum)
%  est_code= bp2(rec, H, sigma,IterNum)  
%    est_code=BP4(rec,H,sigma, IterNum);
%   est_code =BF2(n,k,rec,H,d)
%   est_code= LDPC_decoder(rec,sigma,H,IterNum);  %%%%%%%概率域BP
        est_code = BF(rec1,H,IterNum);
% code=BF(rec,H,IterNum)
%est_code =code(1:k);
        err2=length(find(est_code~=c));
       
        Err2=Err2+err2;
    end
    BER2(nEN2)=Err2/Npf;
    if BER2(nEN2)<1/Npf;
        BER2(nEN2)=.1/Npf;
    end  
  end
%   for nEN3=1:length(EbN0db)
%     en1=10^(EbN0db(nEN3)/10);
%     sigma=1/sqrt(2*rate*en1);
%     nframe1=0;
%     Err3=0;
%     while nframe1<ferrlim
%        nframe1=nframe1+1;
% %        code=c;
% load msg
% 
% 
% %        I=1-2*msg;
%        rec=msg+sigma*randn(1,k);        
% %        [est_code1,success]= LDPC_decoder(rec1,sigma,H,IterNum);  
% %        est_code22=BP4(rec1,H,sigma, IterNum);%%%%%%%概率域BP
%  for i=1:k
%             if rec(i) < 0.5
%          rec(i) = 0;
%       else
%          rec(i) = 1;
%             end
%         end
%     
%        est_code3= rec(1:k);
%        err3=length(find(est_code3~=msg));
%        Err3=Err3+err3;
%     end
%       BER3(nEN3)=Err3/(ferrlim*k);
%      if BER3(nEN3)<1/(ferrlim*k);
%       BER3(nEN3)=.1/(ferrlim*k);
%      end  
%  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------描绘误比特率曲线------------------------------------
figure(3);  
semilogy(EbN0db,BER1,'r-o'); hold on; 
semilogy(EbN0db,BER2,'g-*');   hold on;   %没有采用纠错编码时的BER曲线  hold on;  
% semilogy(EbN0db,BER3,'k-*'); hold on; 
axis([0 Ndb 1/Npf 1])
title('误比特率性能比较');
xlabel('信噪比(dB)');
ylabel('误码率');
 
legend('BP译码','BF译码');
% axis([1 EbN0db a 1]);