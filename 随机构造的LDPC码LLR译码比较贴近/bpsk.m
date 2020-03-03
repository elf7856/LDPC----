function [waveform]=bpsk(bitseq)
for i=1:length(bitseq)
   if bitseq(i)==1
      waveform(i)=-1;
   else
      waveform(i)=1;
   end
end