
clear all;
close all;
sysCfg=sysCfgStr();
global puschDMRS;
global FFTLxL;
subcarriers=sysCfg.subcarriers;
Lengthdelay=sysCfg.maxGroupDelay;

Wfft=exp(-j*2*pi/sysCfg.fftsize);
%PP=2*eNB.Nrb;%RS point number!!!
for p=1:subcarriers/2+1
    FFTpl(p,:)=power(exp(-j*2*pi/sysCfg.fftsize*(p-1)),[0:1:Lengthdelay-1]);
end
for p=subcarriers/2+2:subcarriers+1
    FFTpl(p,:)=power(exp(-j*2*pi/sysCfg.fftsize*(p-1)),[0:1:Lengthdelay-1]);
end
%sys.FFTpl=FFTpl;
FFTLxL=FFTpl'*FFTpl; % ctranspose - Complex conjugate transpose   


%% generate DMRS
ue.NPUSCHID = 42;
ue.NDMRSID = 1;
ue.NSubframe = 0;
chs.PRBSet = (0:24).';
puschDMRS = ltePUSCHDRS(ue,chs);
DMRS=[puschDMRS(1:300).';puschDMRS(301:end).']; % only use the DMRS(1,:)

Ber=[];
%load('Bers.mat');
[txWaveForm,txDmrs,txData]=TxSimu();
for i=1:6
    for SNR=1:13
        BerLS01(SNR,i)=0;
        BerLS02(SNR,i)=0;
        BerLS03(SNR,i)=0;
        for j=1:6
            txWaveFormWithCh=channelpath(txWaveForm,'EPA',SNR);
            %BerRaw(SNR,i)=runonce('noCHE','noEQ',SNR,'awgn');
            %Ber(SNR,i)=runonce('noCHE','noEQ',SNR,'awgn');
            %Ber1(SNR,i)=runonce('noCHE','noEQ',SNR,'rayleigh');

            %BerLS0(SNR,i)=runonce('LS0','ZF',SNR,'awgn');
            %BerLS01(SNR,i)=runonce('LS0','ZF',SNR,'multipath');
            BerLS01(SNR,i)=BerLS01(SNR,i)+ RxOnce(txWaveFormWithCh,txDmrs,txData,'LS0','ZF');
            BerLS02(SNR,i)=BerLS02(SNR,i)+ RxOnce(txWaveFormWithCh,txDmrs,txData,'LS1','ZF');

            %BerLS03(SNR,i)=runonce('LS0','LMMSE',SNR,'raylei');

            %[txWaveFormWithCh,txDmrs,txData]=TxSimu(SNR,'ETU');
            BerLS03(SNR,i)=BerLS03(SNR,i)+ RxOnce(txWaveFormWithCh,txDmrs,txData,'LMMSE','ZF');

            %BerLS1(SNR,i)=runonce('LS1','ZF',SNR,'awgn');
            %BerLS21(SNR,i)=runonce('LS1','ZF',SNR,'multipath');

            %BerLS2(SNR,i)=runonce('LMMSE','ZF',SNR,'awgn');
            %BerLS21(SNR,i)=runonce('LMMSE','ZF',SNR,'multipath');
            %BerLS02(SNR,i)=runonce('LMMSE','ZF',SNR,'multipath');
            %BerLS03(SNR,i)=RxOnce(txWaveFormWithCh,txDmrs,txData,'LMMSE','LMMSE');
        end
        BerLS01(SNR,i)= BerLS01(SNR,i)/j;
        BerLS02(SNR,i)= BerLS02(SNR,i)/j;
        BerLS03(SNR,i)= BerLS03(SNR,i)/j;
    end
    
end

%save Bers.mat Ber Ber1 BerLS0 BerLS01 BerLS1 BerLS11 BerLS2 BerLS21;
%drawResult();
%b=sum(BerRaw,2)/10;
%semilogy(1:SNR,b,'b--');hold on;
maxV=size(BerLS02,1);
dataLen=size(BerLS02,2);

%b=sum(BerRaw,2)/dataLen;
%semilogy(1:maxV,b,'b--');hold on;
% 
% b=sum(BerLS2,2)/dataLen;
% semilogy(1:maxV,b,'g--');hold on;
% 
% b=sum(BerLS01,2)/dataLen;
% semilogy(1:maxV,b,'b-o');hold on;
b=sum(BerLS01,2)/dataLen;
semilogy(1:maxV,b,'b-*');hold on;

b=sum(BerLS02,2)/dataLen;
semilogy(1:maxV,b,'g--');hold on;

b=sum(BerLS03,2)/dataLen;
semilogy(1:maxV,b,'r-o');hold on;
% 
% b=sum(BerLS21,2)/dataLen;
% semilogy(1:maxV,b,'r-o');hold on;

% b=sum(BerLS22,2)/dataLen;
% semilogy(1:maxV,b,'r-*');hold on;

semilogy(maxV:maxV,0.01:0.1:0.01,'g.');grid on;


%save Bers.mat Ber Ber1 BerLS0 BerLS01 BerLS1 BerLS11 BerLS2 BerLS21;

%scatterplot(RxDataTd(1,:))

%a=txWaveFormWithCh(37:end);% remove CP
%a=fft(a,512);% FFT
%a=ifft(a,300);% IDFT
%scatterplot(a)
%


%'OK'