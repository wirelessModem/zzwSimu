function ber=runonce(CHE,EQ,SNR)
global puschDMRS;
sysCfg=sysCfgStr();
symN=sysCfg.subcarriers;
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
        txSymbols=lteSymbolModulate(randi(2,1,sysCfg.modbits*symN)-1, sysCfg.modm);
        txData=[txData,txSymbols];
    end
    %% do DFT :sysCfg.subcarriers
    %% close DFT
    txPrecodingSym=txSymbols;%fft(txSymbols,sysCfg.subcarriers);
    %if(i==4)
    %    txDmrsFd = txPrecodingSym;
    %end
    %% do IFFT :sysCfg.fftsize [shift processing!!!!]
    a=ifft(txPrecodingSym,sysCfg.fftsize);
    if(i==1)
        WaveForm=[a(end-sysCfg.firstCp+1:end);a];% add CP
    else
        WaveForm=[a(end-sysCfg.normalCp+1:end);a];% add CP
    end
    
    txWaveForm=[txWaveForm;WaveForm];
    
   

end

%% wireless channel
a=fft(txWaveForm); %F2T
txWaveFormWithCh=awgn(a,SNR);%,'measured');
txWaveFormWithCh=ifft(txWaveFormWithCh); % T2F
%txWaveFormWithCh=txWaveForm;

%% Rx processing

%% Data Extract
%Atc 
TimeOffset=0;
a=txWaveFormWithCh;
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

rxDmrsFd=rxDmrs;%fft(rxDmrs,sysCfg.fftsize);
rxDmrsFd=rxDmrsFd(1:sysCfg.subcarriers);% need demapping!!!
for i=1:6
    a=rxData(:,i);%fft(rxData(:,i),sysCfg.fftsize);
    a=a(1:sysCfg.subcarriers);% need demapping!!!
    rxDataFd=[rxDataFd,a];
    %scatterplot(a)
end

%% %%%%%%%%%%%%%channel estimation
%% %%%%%%%%%%%%%CHE LS mothod
if strcmp(CHE,'LS')%CHE=='LS'
    txDmrsFd=ifft(txDmrs,sysCfg.fftsize);
    Hls=rxDmrs./txDmrsFd;

    %Hls=conj(txDmrs).*rxDmrsFd;
    %% time field select
    ht=fft(Hls,sysCfg.subcarriers);
    ht(16:end)=0;
    %ht(abs(ht)<max(abs(ht))/20)=0;
    Hd=fft(ht,sysCfg.subcarriers);
    Hls=Hd;
    %%
    pw=abs(Hls);
    Hd=Hls./pw;
else
    Hd=ones(sysCfg.subcarriers,1);
end

%% %%%%%%%%%%%%%EQ
rxDataTd = [];

if strcmp(EQ,'ZF')
    %% %%%%%%%%%%%%%ZF
    for i=1:6
        aFd = Hd.*rxDataFd(:,i);
        aTd=fft(aFd,sysCfg.subcarriers);
        rxDataTd = [rxDataTd,aTd];
    end
else
    %% %%%%%%%%%%%%%draw raw RxData
    for i=1:6
        aFd = Hd.*rxDataFd(:,i);
        aTd=fft(aFd,sysCfg.subcarriers);
        rxDataTd = [rxDataTd,aTd];
        %scatterplot(aFd)
    end
end


%subplot(2,1,1);
%scatter(real(rxDataTd(:,1)),imag(rxDataTd(:,1)));
%scatterplot(RxDataTd(1,:))

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