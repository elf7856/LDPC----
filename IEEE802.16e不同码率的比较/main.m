close all
clear all
%基本参数
% n=2304;                                                                   %%%总码长
% k=1152;                                                                    %%%信息位长度
n=576;                                                                   %%%总码长
k=288;
k1=384;
k2=432;
BER=0;
BER1=0;
BER2=0;
rate=k/n;%%%码率
rate1=k1/n;
rate2=k2/n;
IterNum=30;
%定义仿真参数
ferrlim=30;
Ndb=6;
EbN0db=0:0.5:Ndb;
Npf=ferrlim*k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for nEN1=1:length(EbN0db)
    en1=10^(EbN0db(nEN1)/10);
    sigma=1/sqrt(2*rate*en1);
    nframe1=0;
    Err=0;
    while nframe1<ferrlim
       nframe1=nframe1+1;
       msg = round(rand(1,k));
       [H,c]=ldpc_matrix(msg);
       code=c;
       I=1-2*code;
       rec=I+sigma*randn(1,n); 
     est_code=BP4(rec,H,sigma, IterNum)
     est_code0= est_code(1:k);
     err=length(find(est_code0~=msg));
      Err=Err+err;
     end
      BER(nEN1)=Err/(ferrlim*k);
     if BER(nEN1)<1/(ferrlim*k);
      BER(nEN1)=.1/(ferrlim*k);
     end  
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for nEN1=1:length(EbN0db)
    en1=10^(EbN0db(nEN1)/10);
 sigma1=1/sqrt(2*rate1*en1);
nframe1=0;
     Err1=0;
    while nframe1 < ferrlim
       nframe1=nframe1+1;
       msg1 = round(rand(1,k1));
       [H1,c1]=ldpc_matrix1(msg1);
       code1=c1;
      I1=1-2*code1;
       rec1=I1+sigma1*randn(1,n);
       est_code11=BP4(rec1,H1,sigma1, IterNum);%%%%%%%概率域BP
       est_code1= est_code11(1:k1);
       err1=length(find(est_code1~=msg1));
      Err1=Err1+err1;
    end
      BER1(nEN1)=Err1/(ferrlim*k1);
     if BER1(nEN1)<1/(ferrlim*k1);
      BER1(nEN1)=.1/(ferrlim*k1);
     end  
 end
 for nEN1=1:length(EbN0db)
    en1=10^(EbN0db(nEN1)/10);
    sigma2=1/sqrt(2*rate2*en1);
    nframe1=0;
    Err2=0;
    while nframe1<ferrlim
       nframe1=nframe1+1;
       msg2 = round(rand(1,k2));
       [H2,c2]=ldpc_matrix2(msg2);
       code2=c2;
       I2=1-2*code2;
       rec2=I2+sigma2*randn(1,n);
       est_code22=BP4(rec2,H2,sigma2, IterNum)
       est_code2= est_code22(1:k2);
       err2=length(find(est_code2~=msg2));
       Err2=Err2+err2;
     end
     BER2(nEN1)=Err2/(ferrlim*k2);
     if BER2(nEN1)<1/(ferrlim*k2);
        BER2(nEN1)=.1/(ferrlim*k2);
     end  
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------描绘误比特率曲线------------------------------------
figure(2);
semilogy(EbN0db,BER,'ro-');   hold on;                        %%%最小和  %译码后的BER曲线
semilogy(EbN0db,BER1,'b+-');   hold on;      %%%概率
semilogy(EbN0db,BER2,'g-*');   hold on;   
axis([0 Ndb 1/Npf 1])
% xlim([0 6]);
title('误比特率性能比较');
xlabel('Eb/N0(dB)');
ylabel('误比特率');
legend('rate=1/2','rate=2/3B','rate=3/4A');%,'对数误码率'

