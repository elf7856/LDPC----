tic
m=9;
n=12;
k=n-m;
R=(n-m)/n;
frame_num=10;
Npf=k*frame_num;
EbN0db=0:0.5:6;
for nEN=1:length(EbN0db)
   Err=0; 
   sigma_2=1/(2*10^(EbN0db(nEN)/10)*R);
    for num=1:frame_num
        num;
        s=round(rand(1,n-m));                 %随机产生长为(n-m)的信息序列
        load G
        %getG();
        c=mod(s*G,2);                         %LDPC编码

        waveform=bpsk(c);                     %BPSK调制           

        y=waveform+sqrt(sigma_2)*randn(1,n);  %加性高斯白噪声信道

        maxiter=100;                        %设置最大译码迭代次数maxiter
        [v]=BP1(y,H,sigma_2,maxiter);      %LDPC译码(概率域(SPA1)和对数域上(SPA2)的和积算法,最小和算法(MSA))

        v0=v(m+1:n);
        Err=Err+length(find(s~=v0));               %寻找错误信息位

    end
                %计算比特误码率BER
    BER(nEN)= Err/Npf
    if BER(nEN)<1/Npf;
       BER(nEN)=.1/Npf;
    end
end
figure(1);
    semilogy(EbN0db,BER,'go-');
    xlabel('信噪比/dB')
    ylabel('误码率')
    grid on
    title('average BER')
    axis([0 6 1/Npf 1])
    legend('LDPC编码')
toc
t=toc

















