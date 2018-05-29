
function [channelOut]=channelpath(Timesignal)
global sys;
channelchoise =sys.Ch;%'awgn';
%channel type {'PedB', 'PedBcorr'}
%ChanMod.PDP_dB = [0   -0.9  -4.9  -8    -7.8  -23.9; % Average power [dB]
%                  0 200*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9]; % delay (s)

%ChanMod.PDP_dB = [0   -0.9  -4.9  -8    -7.8  -23.9; % Average power [dB]
                  %%%0 1*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9]; % delay (s)
%                  0 1*10^-9 80*10^-9 120*10^-9 230*10^-9 370*10^-9]; % delay (s)

ChanMod.PDP_dB = [0   -0.9; % Average power [dB]
                  %%%0 1*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9]; % delay (s)
                  0 2000*10^-9]; % delay (s)
ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));

ChanMod2.PDP_dB = [0; % Average power [dB]
                  %%%0 1*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9]; % delay (s)
                  0]; % delay (s)
ChanMod2.normH = sqrt(sum(10.^(ChanMod2.PDP_dB(1,:)/10)));

switch channelchoise
    case 'awgn'
        channelmodel=ChanMod2;
    case 'multipath1'
        channelmodel=ChanMod;
    case 'raylei'
%         ts=1/30720;
%         x=rayleighchan(ts,1);
%         x.StoreHistory=1;
%         x.PathDelay=[0 10*ts 20*ts];
%         x.AvgPathGaindB=[1 0.8 0];
        
        channelOut0=filter(sys.rayleih1,Timesignal(1,:));
        channelOut1=filter(sys.rayleih2,Timesignal(2,:));
end

if isequal(channelchoise,'raylei') == 0
    %bandwidth=20M
    fader=round(channelmodel.PDP_dB(2,:)*1000*sys.Nfft*15);

    unique_fader=unique(fader);
    %G=[1+randn(1,length(unique_fader))/20+j*(1+randn(1,length(unique_fader))/20);1+randn(1,length(unique_fader))/20+j*(1+randn(1,length(unique_fader))/20)];
    G=[ones(1,length(unique_fader))+j*(ones(1,length(unique_fader)));ones(1,length(unique_fader))+j*(ones(1,length(unique_fader)))]/sqrt(2);
    h=[];
    for fadeidx=1:length(unique_fader)
        curr_fader=(fader==unique_fader(fadeidx));
        h(:,unique_fader(fadeidx)+1)=[sqrt(sum(10.^(channelmodel.PDP_dB(1,curr_fader)./10))).*G(:,fadeidx)];
    end
    h=h./channelmodel.normH;
    %y=[];
    channelOut0 = conv(h(1,:), Timesignal(1,:));

    channelOut1 = conv(h(2,:), Timesignal(2,:));
end


channelOut = [channelOut0(1:sum(sys.Ncplength(1:sys.Osn))+sys.Nfft*sys.Osn);channelOut1(1:sum(sys.Ncplength(1:sys.Osn))+sys.Nfft*sys.Osn)];

if sys.debugflag~=0
    save channelpath;
end


%% Add Noise
n = 10^(-sys.snr/20)*(randn(1,size(channelOut0,2)) + 1i*randn(1,size(channelOut0,2)))/sqrt(2);%;/sqrt(ChanMod.nRX);
channelOut0 = channelOut0 + n*sqrt(2048/1200);% not right for 3M and 1.4M
sys.Channelnoise(1,:)=n;
%channelOut0 = channelOut0 + n*sqrt(1200/2048);
n = 10^(-sys.snr/20)*(randn(1,size(channelOut1,2)) + 1i*randn(1,size(channelOut1,2)))/sqrt(2);%;/sqrt(ChanMod.nRX);
channelOut1 = channelOut1 + n*sqrt(2048/1200);% not right for 3M and 1.4M
sys.Channelnoise(2,:)=n;
%channelOut1 = channelOut1 + n*sqrt(1200/2048);
channelOut = [channelOut0(1:sum(sys.Ncplength(1:sys.Osn))+sys.Nfft*sys.Osn);channelOut1(1:sum(sys.Ncplength(1:sys.Osn))+sys.Nfft*sys.Osn)];

%sys.Channelht=h;






