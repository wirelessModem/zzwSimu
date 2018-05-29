function ber=runonce(CHE,EQ,SNR,chan)
global puschDMRS;
sysCfg=sysCfgStr();
subcarriers=sysCfg.subcarriers;
ts=sysCfg.ts;
fd=0;
fftsize=sysCfg.fftsize;
maxPathNum=sysCfg.maxPathNum;
%% ready Tx data
txData=[];
txWaveForm=[];
for i=1:7
    %% generate Tx modulate symbols
    if(i==4)
        %txDmrs=DMRS(1,:).';
        txDmrs=puschDMRS(1:sysCfg.subcarriers);
        txSymbols=txDmrs;
    else
        txSymbols=lteSymbolModulate(randi(2,1,sysCfg.modbits*subcarriers)-1, sysCfg.modm);
        txData=[txData,txSymbols];
    end
    %% do DFT :sysCfg.subcarriers
    %% close DFT
    txPrecodingSym=txSymbols;%fft(txSymbols,sysCfg.subcarriers);
    %if(i==4)
    %    txDmrsFd = txPrecodingSym;
    %end
    %% Data->SubCarrier Mapping[Shift]
    fftInData=subMapFreq(txPrecodingSym,subcarriers,fftsize);

    %% do IFFT :sysCfg.fftsize [shift processing!!!!]
    a=ifft(fftInData,sysCfg.fftsize);
    if(i==1)
        WaveForm=[a(end-sysCfg.firstCp+1:end);a];% add CP
    else
        WaveForm=[a(end-sysCfg.normalCp+1:end);a];% add CP
    end
    
    txWaveForm=[txWaveForm;WaveForm];
    
   

end

%% wireless channel
%a=fft(txWaveForm); %F2T
%txWaveFormWithCh=awgn(a,SNR);
txWaveFormWithCh=channelpath(txWaveForm,chan,SNR);

% if strcmp(chan,'awgn')
%     %txWaveFormWithCh=a;
%     txWaveFormWithCh=awgn(a,SNR);%,'measured');
% elseif strcmp(chan,'rayleigh') % and +awgn
%     %txWaveFormWithCh=awgn(a,SNR);%,'measured');
%     %txWaveFormWithCh=a;
%     %c1 = rayleighchan(1e-5,130) % Create object.
%     %c1.PathDelays = [0 1e-6]    % Change the number of delays.
%     tau=[0 12*ts];% 5*ts];
%     chan=rayleighchan(ts,fd,tau);  
%     %c1.PathDelays = [0 1e-6] 
%     %chan.PathDelays=[ 0 4*ts];
%     txWaveFormWithCh=filter(chan,a);%txWaveFormWithCh);%¹ýÐÅµÀ  
% end

%txWaveFormWithCh=ifft(txWaveFormWithCh); % T2F
%txWaveFormWithCh=txWaveForm;

%% Rx processing

%% Data Extract
%Atc 
TimeOffset=0;
a=reshape(txWaveFormWithCh,length(txWaveFormWithCh),1);
rxData=[];
for i=1:7
    if(i==1)
        rxData=[rxData,a(sysCfg.firstCp+1+TimeOffset:sysCfg.firstCp+sysCfg.fftsize+TimeOffset)];
        a=a(sysCfg.fftsize+sysCfg.firstCp+1:end);
    elseif(i==4)
        rxDmrs=a(sysCfg.normalCp+1+TimeOffset:sysCfg.normalCp+sysCfg.fftsize+TimeOffset);
        a=a(sysCfg.fftsize+sysCfg.normalCp+1:end);
    else
        rxData=[rxData,a(sysCfg.normalCp+1+TimeOffset:sysCfg.normalCp+sysCfg.fftsize+TimeOffset)];
        a=a(sysCfg.fftsize+sysCfg.normalCp+1:end);
    end
    
end

%% T2F
rxDataFd=[];

rxDmrsFd=fft(rxDmrs,sysCfg.fftsize);
%rxDmrsFd=rxDmrsFd(1:sysCfg.subcarriers);% need demapping!!!
for i=1:6
    %a=rxData(:,i);%fft(rxData(:,i),sysCfg.fftsize);
    a=fft(rxData(:,i),sysCfg.fftsize);
    %a=a(1:sysCfg.subcarriers);% need demapping!!!
    rxDataFd=[rxDataFd,a];
    %scatterplot(a)
end

%% SubCarrier->Data DeMapping[Shift]
rxDmrsSym=FreqMapSub(rxDmrsFd,subcarriers);
%rxDmrsSym=[rxDmrsFd(end-(subcarriers/2)+1:end);rxDmrsFd(2:subcarriers/2+1)];
for i=1:6
    rxDataSym(:,i)=FreqMapSub(rxDataFd(:,i),subcarriers);
    %rxDataSym(:,i)=[rxDataFd(end-(subcarriers/2)+1:end,i);rxDataFd(2:subcarriers/2+1,i)];
end
%scatterplot(rxDataSym(:,1));
    
%% %%%%%%%%%%%%%channel estimation
%% %%%%%%%%%%%%%CHE LS mothod
if strcmp(CHE,'LS0')%CHE=='LS'
    %txDmrsFd=ifft(txDmrs,sysCfg.fftsize);
    %rxDmrsp=fft(rxDmrs,sysCfg.fftsize);
    %Hls=rxDmrsFd./txDmrsFd;
    %Hd=conj(txDmrs).*rxDmrsFd;
    Hls=rxDmrsSym./txDmrs;
    
    %% time field select
    a=subMapFreq(Hls,subcarriers,fftsize);
    h_td=ifft(a);
    %b(32:end)=0;
    h_td(48+1:end-48)=0;
    %b(sysCfg.firstCp+1:end)=0;
    %b(abs(b)<max(abs(b))/10)=0;
    c=fft(h_td);
    Hd=FreqMapSub(c,subcarriers);
    Hd_abs=abs(Hd).^2;
    Hd=conj(Hd)./Hd_abs;
    %ht=ifft(Hls,sysCfg.subcarriers);
    %ht(16:end)=0;
    %ht(abs(ht)<max(abs(ht))/20)=0;
    %Hd=fft(ht,sysCfg.subcarriers);
    %Hls=Hd;
    %%
    %pw=abs(Hls);
    %Hd=Hls./pw;
elseif strcmp(CHE,'LS1')%CHE=='LS'
    %txDmrsFd=ifft(txDmrs,sysCfg.fftsize);
    %rxDmrsp=fft(rxDmrs,sysCfg.fftsize);
    %Hls=rxDmrsFd./txDmrsFd;
    %Hd=conj(txDmrs).*rxDmrsFd;
    %Hls=conj(txDmrs).*rxDmrsSym;
    Hls=rxDmrsSym./txDmrs;
    %% time field select
    a=subMapFreq(Hls,subcarriers,fftsize);
    b=ifft(a);
    %b(32:end)=0;
    windth=sysCfg.firstCp;
    b(windth+1:end-windth)=0;
    %b(abs(b)<max(abs(b))/10)=0;
    c=fft(b);
    Hd=FreqMapSub(c,subcarriers);
    Hd_abs=abs(Hd).^2;
    Hd=conj(Hd)./Hd_abs;
    %Hd=(Hd)./Hd_abs;
    
    %ht=ifft(Hls,sysCfg.subcarriers);
    %ht(16:end)=0;
    %ht(abs(ht)<max(abs(ht))/20)=0;
    %Hd=fft(ht,sysCfg.subcarriers);
    %Hls=Hd;
    %%
    %pw=abs(Hls);
    %Hd=Hls./pw;    
elseif strcmp(CHE,'LMMSE')%CHE=='LS'
    %txDmrsFd=ifft(txDmrs,sysCfg.fftsize);
    %rxDmrsp=fft(rxDmrs,sysCfg.fftsize);
    %Hls=rxDmrsFd./txDmrsFd;
    %Hd=conj(txDmrs).*rxDmrsFd;
    %Hls=conj(txDmrs).*rxDmrsSym;
    %Hmmse=CE_lmmse(Yrs,Nrb,RS,Lengthdelay,ppsMaxPathnum,Nfft)
    
    Hls=rxDmrsSym./txDmrs;
    
    %% time field select
    a=subMapFreq(Hls,subcarriers,fftsize);
    h_td=ifft(a);
    %b(32:end)=0;
    h_td(48+1:end-48)=0;
    %b(sysCfg.firstCp+1:end)=0;
    %b(abs(b)<max(abs(b))/10)=0;
    c=fft(h_td);
    Hd=FreqMapSub(c,subcarriers);
    Hd_abs=abs(Hd).^2;
    Hd1=conj(Hd)./Hd_abs;
    %figure;plot(real(Hd1),'b--');hold on;
    
    
    Hd=CE_lmmse(rxDmrsSym,sysCfg.Nrb,txDmrs,sysCfg.maxGroupDelay,sysCfg.maxPathNum,sysCfg.fftsize);
    %Hd=conj(Hd)/sqrt(2);
    Hd_abs=abs(Hd).^2;
    %Hd_abs=Hd_abs*6;
    %Hd=Hd./Hd_abs;%
    Hd=conj(Hd)./Hd_abs;% result is same..?? +conj is right
    
%     figure;plot(real(Hd1),'b--');hold on;
%     plot(real(Hd),'r--');
%     figure;plot(imag(Hd1),'b--');hold on;
%     plot(imag(Hd),'r--');
    
elseif strcmp(CHE,'LMMSExx')
    Hlsx=conj(txDmrs).*rxDmrs
else
    Hd=ones(sysCfg.subcarriers,1);
end

%% %%%%%%%%%%%%%EQ
rxDataTd = [];

if strcmp(EQ,'ZF')
    %% %%%%%%%%%%%%%ZF
    for i=1:6
        aFd = Hd.*rxDataSym(:,i);
        %aFd = conj(Hd).*rxDataSym(:,i);
        aTd=aFd;%fft(aFd,sysCfg.subcarriers);
        rxDataTd = [rxDataTd,aTd];
    end
else
    %% %%%%%%%%%%%%%draw raw RxData
    for i=1:6
        aFd = Hd.*rxDataSym(:,i);
        aTd=aFd;%fft(aFd,sysCfg.subcarriers);
        rxDataTd = [rxDataTd,aTd];
        %scatterplot(aFd)
    end
end

%scatterplot(rxDataTd(:,1));
%subplot(2,1,1);
%scatter(real(rxDataTd(:,1)),imag(rxDataTd(:,1)));
% a=Hd1.*rxDataSym(:,1);
% scatter(real(a(:,1)),imag(a(:,1)));hold on;
% scatter(real(rxDataTd(:,1)),imag(rxDataTd(:,1)));
% close all;
%% demod
a=[];
for i=1:2^sysCfg.modbits
    a=[a;uint8(dec2binvec(i-1,sysCfg.modbits)')];
end
modTbl=lteSymbolModulate(a, sysCfg.modm);

rxData=[];
for i=1:size(rxDataTd,2)
    for j=1:size(rxDataTd,1)
        [x,y]=min(rxDataTd(j,i)-modTbl);
        rxData(j,i)=modTbl(y);
        if rxData(j,i)==0
            rxData(j,i)=rxData(j,i);
        end
    end
    
end
 

%rxDataTd = [];
a=sum(sum(abs(txData-rxData)~=0));
ber=a/(size(rxData,1)*size(rxData,2));
%subplot(2,1,2);scatter(real(rxDataTd(:,1)),imag(rxDataTd(:,1)));
end