clear all
close all
clc
%-----------------------------升余弦滤波器
Fs=9600;            %采样频率
Ts=1/Fs;            %采样间隔
Fd=960;              %Doppler频偏，以Hz为单位
tau=[0,0.002];          %多径延时，以s为单位
pdf=[0,0];          %各径功率，以dB位单位
h=rayleighchan(Ts,Fd,tau,pdf);
%-------------------------------通过信道
data1=[1,zeros(1,100)];%数据1 
fc=96;
t=0:1/Fs:999*(1/Fs);
data2=cos(2*pi*fc*t); %数据2 
data2_fft=fft(data2,100);
data2_abs=abs(data2_fft);

y1=filter(h,data1);
y11=abs(y1);
y2=filter(h,data2);
y22=fft(y2,100);
y222=abs(y22);

subplot(2,2,1);
x1=1:length(data1);
h1=stem(x1,data1);
set(h1,'MarkerFaceColor','red')

subplot(2,2,3);
x3=1:length(y11);
h3=stem(x3,y11);
set(h3,'MarkerFaceColor','red')

subplot(2,2,2);
x2=1:length(data2_abs);
h2=stem(x2,data2_abs);
set(h2,'MarkerFaceColor','blue')

subplot(2,2,4);
x4=1:length(y222);
h4=stem(x4,y222);
set(h4,'MarkerFaceColor','blue')