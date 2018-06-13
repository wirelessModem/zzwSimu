classdef sysCfgStr
    properties(Constant=true)
        %subcarrier 15KHz
        subcarriers=300;%%300;
        bw=5e6;
        fftsize=512;
        firstCp=40;
        normalCp=36;
        maxPathNum=12;
        maxGroupDelay=40;%floor(512/6);%44;%36;%36;%150;
        samplerate=7.68e6;
        ts=1/7.68e6;
        tsNsec=(10^9)*(1/7.68e6);
        symboln=7.68e6/512/1000; %ms
        Nrb=25;
        %fftsize=256;
        %firstCp=20;
        %normalCp=18;
        %samplerate=3.84e6;
        %ts=1/3.84e6;
        %symboln=3.84e6/256/1000; %ms
        
        
        modm='16QAM';%'QPSK';%'16QAM';%[1:BPSK,2:QPSK,4:16QAM,6:64QAM]
        modbits=4;%2;%4;
        
    end
end