classdef sysCfgStr
    properties(Constant=true)
        %subcarrier 15KHz
        bw=5e6;
        fftsize=512;
        samplerate=7.68e6;
        ts=1/7.68e6;
        symboln=7.68e6/512/1000; %ms
    end
end