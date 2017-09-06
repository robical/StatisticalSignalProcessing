%% Deconvoluzione MMSE nelle frequenze con filtro di Wiener

clear all;
close all;
clc;

%Parametri

N=201;
sigw=10;
sigz=1;
e=1e-2; %quanto è vicino al cerchio unitario lo zero?
Nmed=500; %numero di iterazioni Monte Carlo

puls_norm=-pi+pi/(N/2):pi/(N/2):pi; %asse della pulsazione normalizzata

Rx=zeros(N,1); %inizializzazione della trasformata Z dell'autocorr
Rw=zeros(N,1); %inizializzazione della trasformata Z dell'autocorr
Rz=zeros(N,1);

for kk=1:Nmed

%processo bianco in ingresso con potenza sigw

w=sqrt(sigw)*randn(N,1);

%w=cos(2*pi*10*[1:N]');

W=fftshift(fft(w)); %trasformata centrata in zero
RW=(1/N)*W.*conj(W);
%RW=diag(RW); %trasformata di Fourier dell'autocorrelazione (non mediata)

Rw=Rw+RW;

%disturbo bianco z con potenza sigz

z=sqrt(sigz)*randn(N,1);

Z=fftshift(fft(z)); %trasformata centrata in zero
RZ=(1/N)*Z.*conj(Z);
Rz=Rz+RZ;
%filtro

Bp=[(1-e)*exp(-1i*(pi/4)) (1-e)*exp(1i*(pi/4))];

A=zero(Bp);

%osservazioni rumorose

temp=conv(w,A); %lunga N+2
temp=temp(2:end-1); %tagliato a N

x=temp+z; %osservazioni rumorose

%trasformata dell'autocorrelazione

X(:,kk)=fftshift(fft(x)); %trasformata centrata in zero

RX=(1/N)*X(:,kk).*conj(X(:,kk));
%RX=diag(RX); %trasformata di Fourier dell'autocorrelazione (non mediata)

Rx=Rx+RX;

end

%Devo conoscere il filtro

H=fftshift(fft(A,N)); %trasformata di Fourier del filtro, centrato in zero
num=sigw*fliplr(conj(H)); %numeratore del filtro di Wiener

%Calcolo delle trasformate delle autocorrelazioni
Rx=Rx/Nmed;
Rw=Rw/Nmed;
Rz=Rz/Nmed;

%Calcolo filtro di Wiener in frequenza
Wi=num'./Rx;

%Calcolo filtro inverso in frequenza

Hinv=1./H;

%Segnale ricevuto deconvoluto - trasformata
% W_stim = con filtro Wiener
% W_inv = on filtro inverso

W_stim=zeros(N,1);
W_inv=zeros(N,1);

for kk=1:Nmed
    W_stim=W_stim+Wi.*X(:,kk);
    W_inv=W_inv+Hinv'.*X(:,kk);
end

W_stim=W_stim/Nmed;
W_inv=W_inv/Nmed;

%errore tra spettro di potenza dell'ingresso "vero" e spettro stimato con
%deconvoluzione Wiener

err_spectr= abs(abs(W_stim).^2-Rw);

MSE=sum(err_spectr); %MSE in frequenza

%Proviamo la stima dello spettro di potenza del segnale in ingresso


% figure,...
%     subplot(1,3,1),plot(puls_norm,Rx,'r'),xlabel('Pulsazione normalizzata [rad]'),...
%     ylabel('Potenza'),title('Densità spettrale di potenza del processo colorato'),...
%     subplot(1,3,2),plot(puls_norm,Rw,'k'),xlabel('Pulsazione normalizzata [rad]'),...
%     ylabel('Potenza'),title('Densità spettrale di potenza del processo ingresso'),...
%     subplot(1,3,3),plot(puls_norm,abs(W_stim).^2,'k'),xlabel('Pulsazione normalizzata [rad]'),...
%     ylabel('Potenza'),title('Densità spettrale di potenza del processo deconvoluto'),...

figure,subplot(1,2,1),plot(puls_norm,abs(W_stim).^2,'r'),xlabel('Pulsazione normalizzata [rad]'),...
    ylabel('Potenza'),title('Ingresso (spettro di potenza) stimato con Wiener'),...
    subplot(1,2,2),plot(puls_norm,abs(W_inv).^2,'r'),xlabel('Pulsazione normalizzata [rad]'),...
    ylabel('Potenza'),title('Ingresso (spettro di potenza) stimato con filtro inverso'),...
    
    
figure,plot(puls_norm,err_spectr,'k'),...
    title(strcat('Densità spettrale di potenza dell''errore della stima per SNR pari a ',num2str(10*log10(sigw/sigz)),'dB')),...
    xlabel('Pulsazione normalizzata [rad]'),ylabel('Potenza'),

disp(strcat('MSE della stima :',num2str(MSE)));
% figure,plot(puls_norm,abs(H),'g'),xlabel('Pulsazione normalizzata [rad]'),...
%     ylabel('Potenza'),title('Funzione di trasferimento in ampiezza del filtro')
%     
    