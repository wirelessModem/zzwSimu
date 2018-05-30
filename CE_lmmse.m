function Hmmse=CE_lmmse(Yrs,Nrb,RS,Lengthdelay,ppsMaxPathnum,Nfft)
global FFTLxL;
subCars = Nrb*12;
%extract data
%Ydata = Ydata(dataPos);

%% get the LS H
Hlsx=(conj(RS).*Yrs);

%% here, we place LS.H to right position without Shift postprocessing!!!
%Hls=subMapFreq(Hlsx,subCars,Nfft);
%k=1:6:subCars;
Hls=[Hlsx(1:subCars/2);0;Hlsx(subCars/2+1:subCars)];

%Hls(6*k-5)=Hlsx(6*k-5);
%Hls=Hlsx;%[Hlsx(1:150);0;0;0;Hlsx(151:end)];

%k=1:6:subCars;
%Hls=[];
%Hls(k)=Hlsx(k);
%% LS.H->time.h
%htime=sqrt(2048)*sqrt(2048/1200)*ifft(Hls,2048);
htime=ifft(Hls,Nfft);
htime=reshape(htime,1,length(htime));
%htime=ifft(Hls,2048);
ht=htime(1:Lengthdelay);%256);%Lengthdelay);
%ht=[htime(end-16+1:end) htime(1:Lengthdelay-16)];

% len=size(Hlsx,1);
% n1=ones(len,1);
% n1=n1*0.000000000000000001i;%Just to ensure that the function awgn adds 'complex gaussian noise'..
% noise=awgn(n1,10);
% variance=var(noise);
% 
% htt=ifft(Hlsx,Nfft);
% httt=htt(1:Lengthdelay);
% a=sort(abs(httt),'descend');
% PthHt=httt(  (abs(httt)>=a(ppsMaxPathnum) ) );%& abs(ht)> (max(abs(ht))/15) );
% Hlss=fft(PthHt,Nfft);
% hh=diag(Hlss(1:subCars));
% hh_myu=sum(hh,1)/subCars;
% hh_mid=hh-hh_myu(ones(subCars,1),:);
% sum_hh_mid=sum(hh_mid,1);
% Rhh=(hh_mid'*hh_mid- (sum_hh_mid'*sum_hh_mid)/subCars)/(subCars - 1);
% Hlmmse=Rhh*inv(Rhh+variance*inv(RS*RS'))*Hlsx;
% Hlmmse=Hlmmse*(Nfft);
% Hlmmse=reshape(Hlmmse,length(Hlmmse),1);

%% normalization!!
%ht=ht./sqrt(sum(ht*ht'));

ppsMaxPathnum=Lengthdelay;
hTmp=sort(abs(ht),'descend');
PthHt=ht(  (abs(ht)>=hTmp(ppsMaxPathnum) ) );%& abs(ht)> (max(abs(ht))/15) );

PathIndex=[1:1:Lengthdelay];
PathIndex=PathIndex( (abs(ht)>=hTmp(ppsMaxPathnum) ) );%& abs(ht)> (max(abs(ht))/15));

PowHt=ht;%(abs(ht)>= (hTmp(ppsMaxPathnum)/2) );
FFTL=FFTLxL(PathIndex,PathIndex);
%FFTL=FFTL/(2*Nrb);
%[V,D]=eig(PthHt'*PthHt);
Rhh=diag(PthHt'*PthHt);
a=ppsMaxPathnum;%size(Rhh,1);
Rhh=diag(Rhh);%diag(ones(a,1));%D,x;%diag(Rhh);%diag(ones(a,1));D;%diag(Rhh)~==diag(ones(a,1));
%Rhh=inv(Rhh);
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
    factor1=8;%8;%10;%15;%25;%96;%25;
    sigma_1=sigma*Nrb*factor1;
    
%else
    %factor2=sqrt(2048/1200)*2; % for awgn is good
    %factor2=sqrt(2048/1200)*3.6; % 
    factor2=4;%8;%4;%sqrt(2048/1200)*4.6; % 
    sigma_2=max(max(Rhh))*factor2;
%    sigma=sys.sigma(2);
sigma=max( sigma_1,sigma_2 );
if (sigma*100)>0.01
    sigma=sigma_1;
else
    sigma=sigma_2;
    
end
%sigma=max( sigma_1,sigma_2 );
deltaDiag=0;
%sigma2=80;
%sigma=sigma_1;
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
Hmmse=[H(1:subCars/2) H(subCars/2+2:subCars+1)];%;%
%Hmmse=FreqMapSub(H,subCars);
Hmmse=Hmmse*(Nfft);
Hmmse=reshape(Hmmse,length(Hmmse),1);