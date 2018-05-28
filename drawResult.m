function drawResult()
load('Bers.mat');
maxV=length(Ber);
dataLen=size(Ber,2);
b=sum(Ber,2)/dataLen;
semilogy(1:maxV,b,'b-o');hold on;
b=sum(Ber1,2)/dataLen;
semilogy(1:maxV,b,'b-*');hold on;

b=sum(BerLS0,2)/dataLen;
semilogy(1:maxV,b,'g-o');hold on;
b=sum(BerLS01,2)/dataLen;
semilogy(1:maxV,b,'g-*');hold on;

b=sum(BerLS1,2)/dataLen;
semilogy(1:maxV,b,'r-s');hold on;
b=sum(BerLS11,2)/dataLen;
semilogy(1:maxV,b,'r-d');hold on;


b=sum(BerLS2,2)/dataLen;
semilogy(1:maxV,b,'r-o');hold on;
b=sum(BerLS21,2)/dataLen;
semilogy(1:maxV,b,'r-*');hold on;

semilogy(maxV:maxV,0.01:0.1:0.01,'g.');grid on;
end