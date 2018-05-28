function Hmmse=CE_lmmse(Yrs,Nrb,RS,Lengthdelay,ppsMaxPathnum,Nfft)
global FFTLxL;
subCars = Nrb*12;
%extract data
%Ydata = Ydata(dataPos);

%% get the LS H
Hlsx=(conj(RS).*Yrs);

%% here, we place LS.H to right position without Shift postprocessing!!!
Hls=Hlsx;%subMapFreq(Hlsx,subCars,Nfft);
%k=1:6:subCars;
%Hls=[Hlsx(1:6:Nrb*6); 0 ;Hlsx(Nrb*6+6:6:end)];

%Hls(6*k-5)=Hlsx(6*k-5);
%Hls=Hlsx;%[Hlsx(1:150);0;0;0;Hlsx(151:end)];

%% LS.H->time.h
%htime=sqrt(2048)*sqrt(2048/1200)*ifft(Hls,2048);
htime=ifft(Hls,Nfft);
htime=reshape(htime,1,length(htime));
%htime=ifft(Hls,2048);
ht=htime(1:Lengthdelay);%256);%Lengthdelay);
%ht=[htime(end-16+1:end) htime(1:Lengthdelay-16)];
%% normalization!!
%ht=ht./sqrt(sum(ht*ht'));


hTmp=sort(abs(ht),'descend');
PthHt=ht(  (abs(ht)>=hTmp(ppsMaxPathnum) ) );%& abs(ht)> (max(abs(ht))/15) );

PathIndex=[1:1:Lengthdelay];
PathIndex=PathIndex( (abs(ht)>=hTmp(ppsMaxPathnum) ) );%& abs(ht)> (max(abs(ht))/15));

PowHt=ht;%(abs(ht)>= (hTmp(ppsMaxPathnum)/2) );
FFTL=FFTLxL(PathIndex,PathIndex);
%FFTL=FFTL/(2*Nrb);

Rhh=diag(PthHt'*PthHt);
Rhh=diag(Rhh);

totalPow=PowHt*PowHt';%sum(abs(PowHt).^2);
if(length(PowHt)>2*ppsMaxPathnum)
    %sigma=totalPow-sum(abs(PthHt).^2)/(length(PowHt)-ppsMaxPathnum);
    sigma=(totalPow-PthHt*PthHt')/(length(PowHt)-length(PthHt));
else
    powTmp=hTmp(1:fix(length(PowHt)/2));
    %sigma=totalPow-sum(abs(powTmp).^2)/fix(length(PowHt)/2);
    sigma=(totalPow-powTmp*powTmp')/fix(length(PowHt)/2);
end

%if sigma>0.0035
    factor1=96;%25;
    sigma_1=sigma*Nrb*factor1;
    
%else
    %factor2=sqrt(2048/1200)*2; % for awgn is good
    %factor2=sqrt(2048/1200)*3.6; % 
    factor2=4;%sqrt(2048/1200)*4.6; % 
    sigma_2=max(max(Rhh))*factor2;
%    sigma=sys.sigma(2);
if (sigma*100)>0.001
    sigma=sigma_1;
else
    sigma=sigma_2;
    
end
%sigma=max( sigma,max(max(Rhh))*factor2 );
deltaDiag=0;
%sigma2=80;
iRR=inv(Rhh);
R=FFTL+(sigma)*iRR+deltaDiag;


%Hin=inv(R)*PthHt.';
iR=inv(R);
for idx=1:size(PthHt,2)
    Hin(idx,1)=(iR(idx,:)*PthHt.');%/abs(sum(iR(idx,:)));
end
Hin=Hin.';
%Hin=Hin/sum(abs(Hin));
hmmse(PathIndex)=Hin;%time Hmmse!!!

%hmmse=[hmmse(17:end),zeros(1,Nfft-length(hmmse)),hmmse(1:16)];

%% fft time.h->Freq.H
H=fft(hmmse,Nfft);%/sqrt(2048/1200);%/sqrt(2048);
%Hmmse=[H(1:subCars/2) H(subCars/2+2:subCars+1)];%[H(1:150),H(152:301)];
Hmmse=H(1:subCars);%FreqMapSub(H,subCars);
Hmmse=Hmmse*(Nfft);
Hmmse=reshape(Hmmse,length(Hmmse),1);