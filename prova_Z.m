%% Implementazione numerica dei poli con equazioni alle differenze
close all;
clear all;
clc;


e=1e-6;

%apro=wavread('prova.wav');



%Sfasatore puro

Bp=[(1-e)*exp(-1i*(pi/4)) (1-e)*exp(1i*(pi/4))]; %sequenza di poli - tutti a fase minima
Az=conj(1./Bp); %sequenza di zeri - di conseguenza a fase max


[x]=zero(Az);
% 
% aud_fil=conv(apro,x);
% 



fc=2*1e3;
fN=fc/2;

[Spettro,f]=spet_plo(Az,Bp,fc);

%effetto su sinusoide campionata

t=0:(1/fc):0.5;

f1=fN/4;
si=sin(2*pi*f1*t);

%Convoluzione con lo zero

si_fil=conv(si,x);

%figure,plot(si_fil),title('Sinusoide a fN/4 filtrata con zero a +-fN/4')

%convoluzione con entrambi i poli (per annullare l'effetto dei 2 zeri)

[y]=polo(si_fil,Bp(1));

[y]=polo(y,Bp(2));

%plotting

%figure,plot(imag(y),'b')

figure,subplot(1,2,1),plot(f,10*log10(abs(Spettro))),title('Modulo della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Ampiezza [dB]'),axis([-fN fN min(10*log10(abs(Spettro)))-1 max(10*log10(abs(Spettro)))+1]),...
    subplot(1,2,2),plot(f,angle(Spettro)),title('Fase della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Fase [rad]'),axis([-fN fN -pi pi]),...


