function [txWaveFormWithCh,txDmrs,txData]=TxSimu(SNR,chan)
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

end