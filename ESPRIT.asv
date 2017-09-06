%% ESPRIT (Estimation of signal parameters via rotational invariance)

clear all;
close all;
clc;

%Il metodo si basa sull'invarianza per traslazione di 1 campione per la
%matrice del sottospazio di segnale, che si ricava come sempre dalla
%decomposizione della matrice di correlazione campionaria delle
%osservazioni

%Ampiezza e fasi variano casualmente

M=3; %sinusoidi complesse in gioco

N=200; %numero di campioni dell'osservazione

sign=2; %Potenza totale del rumore gaussiano circolare complesso
Amp=5:5:25;
SNR=20*log10(Amp/sign);
%gamma_inv=eye(L)*(1/sign); %matrice di covarianza di rumore 



%chiaramente devo definire pulsazioni normalizzate comprese tra -pi e pi
w1=pi/2;
w2=pi/4;
w3=pi*(2/3);

vere=[w1;w2;w3];
%matrice A(w)

A=[exp(-1i*w1*[0:N-1]') exp(-1i*w2*[0:N-1]') exp(-1i*w3*[0:N-1]')];

%osservazioni

for gg=1:length(Amp)
s=raylrnd(repmat(Amp(gg),M,1)).*exp(1i*(2*pi)*rand(M,1)); %vettore ampiezze complesse
%delle sinusoidi

x=(A*s)'+ (sqrt(sign/2))*(randn(1,N)+1i*randn(1,N)); %modello osservazioni con rumore

%Stima della matrice di autocorrelazione campionaria

Rx_s=xcorr(x,'biased'); %1 modo
%Rx_s_trunc=(1/N)*cconv(x,conj(fliplr(x))); %2 modo

N_aut=length(Rx_s); %dispari

rx_semi=Rx_s((N_aut-1)/2+1:end); %solo parte da 0 a N

%Stima campionaria Matrice di covarianza di segnale completa
for i=1:N
    RX(i,:)=Rx_s(N+(i-1):-1:i); %riempio una riga della matrice
end

%Decomposizione autovalori/autovettori
[V,D]=eig(RX);

%Prendo solo quelli afferenti allo spazio di segnale
Us=V(:,1:M); % NxM

%Ora creo le 2 matrici dalle quali cercherò l'invarianza rotazionale
slice1=[eye(N-1) zeros(N-1,1)]; %(N-1)xN
slice2=[zeros(N-1,1) eye(N-1)];

Ub1=slice1*Us; 
Ub2=slice2*Us;

%Trovo la matrice phi, di rotazione, che trasforma Ub1 in Ub2!

phi=(Ub2'*Ub2)\(Ub2'*Ub1); % (M)x(M)

%Ho già gli autovalori, ne prendo la fase e questa è la stima di frequenza
phi_d=diag(phi);

wstim(:,gg)=-angle(phi_d);
end

%figure,plot(SNR,abs(wstim-repmat(vere,1,length(Amp))).^2);

%disp(strcat('Stima delle pulsaz norm con SNR pari a: ',num2str(SNR),' dB'));




%Usare TLS (Total Least square)

