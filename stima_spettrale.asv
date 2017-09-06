%% Prove stima spettrale non parametrica (periodogramma)

clear all;
close all;
clc;

N=300; %numero campioni della realizzazione del processo casuale bianco gaussiano con potenza sigx


sigx=1;

A=1;

xr=sqrt(sigx)*randn(1,N); %l'ampiezza della cosinusoide sarà
%N/2, dove però N è la dimensione del SEMISPETTRO, cioè 400

%filtro FIR ordine 2

z1=0.5;

[Bz]=zero([z1*exp(1i*pi/2) (1/z1)*exp(-1i*pi/2)]);

%Colorato

y=cconv(xr,Bz,N);

aut=(1/N)*cconv(y,conj(fliplr(y)),N);

per2=fftshift(fft(aut));

periodo=fftshift(fft(y));

fk=assefr(N,1); %asse delle freq normalizzate

figure,subplot(1,2,1),plot(fk,(1/N)*abs(periodo).^2,'r'),title('Periodogramma da DFT'),...
    subplot(1,2,2),plot(fk,abs(per2),'g'),title('Periodogramma da stima autocorr circolare')

%Stima MA

%suppongo filtro a supporto limitato 1+Nb, in questo caso Nb=2 perchè lo
%so,divido la sequenza y in M sottosequenze da Nb+1 campioni, delle quali
%calcolo autocorrelazione circolare, e che poi vado a mediare tra loro;
%infine di queste vedo la rappresentazione nel piano Z, che dovrebbe avere
%Nb zeri

Nze=3; %Nb+1
%stimo autocor circ per ogni sottoseq
for i=1:Nze:N-Nze
    sseq(i,:)=y(i:i+2);
    aut_sseq(i,:)=cconv(sseq(i,:),conj(fliplr(sseq(i,:))),Nze);
end

M=size(sseq,1); %numero di sottosequenze estratte

aut_med=(1/M)*sum(aut_sseq,1);


