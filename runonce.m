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
a=fft(txWaveForm); %F2T
if strcmp(chan,'awgn')
    txWaveFormWithCh=awgn(a,SNR);%,'measured');
elseif strcmp(chan,'rayleigh') % and +awgn
    txWaveFormWithCh=awgn(a,SNR);%,'measured');
    %c1 = rayleighchan(1e-5,130) % Create object.
    %c1.PathDelays = [0 1e-6]    % Change the number of delays.
    chan=rayleighchan(ts,fd);  
    txWaveFormWithCh=filter(chan,txWaveFormWithCh);%¹ýÐÅµÀ  
end




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
    windth=sysCfg.firstCp*2;
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
elseif strcmp(CHE,'LS2')%CHE=='LS'
    %txDmrsFd=ifft(txDmrs,sysCfg.fftsize);
    %rxDmrsp=fft(rxDmrs,sysCfg.fftsize);
    %Hls=rxDmrsFd./txDmrsFd;
    %Hd=conj(txDmrs).*rxDmrsFd;
    %Hls=conj(txDmrs).*rxDmrsSym;
    %Hd=CE_lmmse();
    Hls=rxDmrsSym./txDmrs;
    %% time field select
    a=subMapFreq(Hls,subcarriers,fftsize);
    b=ifft(a);
    b(sysCfg.maxGroupDelay+1:end-sysCfg.maxGroupDelay)=0;
    %b=b(1:sysCfg.maxGroupDelay);
    bb=sort(abs(b),'descend');
    b(abs(b)<abs(bb(maxPathNum)))=0;
    %b(32:end)=0;
    %b(abs(b)<max(abs(b))/10)=0;
    %[envHigh, envLow] = envelope(b,4,'peak');
    %envMean = (envHigh+envLow)/2;
    %a=ones(1,4)/4;
    %fb=filter(a,1,b);
    c=fft(b,fftsize);
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
elseif strcmp(CHE,'LMMSE')
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