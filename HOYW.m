%% Metodo HOYW per stima di frequenza

clear all;
close all;
clc;

%Ampiezza e fasi variano casualmente

M=3; %sinusoidi complesse in gioco
N=100; %numero di campioni dell'osservazione
L=10*M; %cioè 9 in questo caso, ordine del predittore HOYW


sign=25:-1:2; %Potenza totale del rumore gaussiano circolare complesso
Amp=25;
SNR=20*log10(Amp./sign);
Mitr=10;

%chiaramente devo definire pulsazioni normalizzate comprese tra -pi e pi
w1=pi/2;
w2=pi/4;
w3=pi*(2/3);

%disp('Pulsazioni normalizzate VERE, oggetto della stima: ');
vere=[w1 w2 w3];

vere=sort(vere);

MaxL=10;

phi1=0.1;
phi2=0.2;
phi3=0.3;



w_svd_tot=zeros(length(sign),M);
w_LS_tot=zeros(length(sign),M);
for itr=1:Mitr
for gg=1:length(sign)

%x=zeros(1,N);    

s=repmat(Amp,M,1).*exp(1i*(2*pi)*[phi1;phi2;phi3]); %vettore ampiezze complesse
%delle sinusoidi    
%matrice A(w)

A=[exp(1i*w1*[0:N-1]') exp(1i*w2*[0:N-1]') exp(1i*w3*[0:N-1]')];
%osservazioni

x= (A*s)'+ (sqrt(sign(gg)/2))*(randn(1,N)+1i*randn(1,N)); %modello osservazioni con rumore

%Stima della matrice di autocorrelazione campionaria

Rx_s=xcorr(x,'biased'); %1 modo
%Rx_s_trunc=(1/N)*cconv(x,conj(fliplr(x))); %2 modo

N_aut=length(Rx_s); %dispari

rx_semi=Rx_s((N_aut-1)/2+1:end); %solo parte da 0 a N

%figure,subplot(1,2,1),plot(abs(Rx_s((N_aut-1)/2+1:end)),'b'),...
%    title(strcat('Stima autocorrelazione campionaria da lag 1 a ',num2str(N)))
    %subplot(1,2,2),plot(abs(Rx_s_trunc),'r')
    
%Devo scegliere ordine del filtro L, dato che il predittore è lungo M


%costruisco la matrice di autocorrelazione
for i=1:L
    RX(i,:)=rx_semi(L+(i-1):-1:i); %riempio una riga della matrice
    rx_v(i)=rx_semi(L+i); %riempio un pezzo del vettore
end

%inverto con decomposizione QR
[Q,R]=qr(RX);
b_v=-Q'*rx_v';
a_stim=R\b_v; %vettore con i coefficienti del filtro, chiaramente va aggiunto l'1 all'inizio

%Inversione con LS
%gamma_inv=eye(size(RX,1))*(1/sign(gg)); %matrice di covarianza di rumore

a_LS=-(RX'*RX)\(RX'*rx_v'); %inversione con inversa generalizzata

%sequenza della quale voglio calcolare i poli (radici)
a_LS=[1;a_LS];

%Calcolo i poli
rad_LS=roots(a_LS);

%Devo scegliere gli M più vicini al cerchio unitario, ovvero quelli con il
%valore del modulo più vicino a 1

%Faccio 3 check, ad ognuno scelgo quello più vicino al cerchio unitario
for i=1:M
    [val,ind]=min(abs(abs(rad_LS)-1));
    w_LS(gg,i)=-angle(rad_LS(ind));
    rad_LS=[rad_LS(1:ind-1);rad_LS(ind+1:end)];
end

%disp('Soluzioni con inversione LS: ');
w_LS(gg,:);
w_LS(gg,:)=sort(w_LS(gg,:)); %dal più piccolo al più grande




%sequenza della quale voglio calcolare i poli (radici)
a_stim=[1;a_stim];

%Calcolo i poli
radix=roots(a_stim);

%Devo scegliere gli M più vicini al cerchio unitario, ovvero quelli con il
%valore del modulo più vicino a 1

%Faccio 3 check, ad ognuno scelgo quello più vicino al cerchio unitario
for i=1:M
    [val,ind]=min(abs(abs(radix)-1));
    wstim(gg,i)=-angle(radix(ind));
    radix=[radix(1:ind-1);radix(ind+1:end)];
end

%disp('Soluzioni con inversione QR: ');
wstim(gg,:);
wstim(gg,:)=sort(wstim(gg,:));

%funziona..

%Altro modo: stima della matrice di autocorrelazione delle osservazioni con
%N>L e decomposizione SVD invertendo con solo gli M autovettori
%corrispondenti agli autovalori maggiori

%Devo scegliere ordine del filtro L, dato che il predittore è lungo M
L=3*M; %cioè 9 in questo caso
Nc=round(N/3); %N>L, problema sovradeterminato
%costruisco la matrice di autocorrelazione
clear RX;
clear rx_v;
for i=1:Nc
    RX(i,:)=x(L+(i-1):-1:i); %riempio una riga della matrice NxL
    rx_v(i)=x(L+i); %riempio un pezzo del vettore Nx1
end

%decomposizione SVD di RX
[U S V]=svd(RX); %V base ortonormale nello spazio LxL
%U base ortonormale nello spazio NxN

Us=U(:,1:M); %stima della base del sottospazio di segnale, è una stima perchè deriva 
%dalla stima campionaria dell'autocorrelazione

Vs=V(:,1:M); %per autovettori sinistri
Lam_s=S(1:M,1:M); %estraggo la parte afferente agli autovalori dello spazio di segnale, che contengono ovviamente
stim_rum=sum(diag(S(M+1:end,M+1:end))); %stima del rumore
Lam_s=Lam_s-stim_rum*eye(size(Lam_s));

%inverto i parametri direttamente con queste matrici
a=-Vs*(Lam_s\Us')*rx_v'; %ecco la sequenza stimata, è gia Lx1

%sequenza della quale voglio calcolare i poli (radici)
a=[1;a];

%Calcolo i poli
rad_svd=roots(a);

%Devo scegliere gli M più vicini al cerchio unitario, ovvero quelli con il
%valore del modulo più vicino a 1

%Faccio 3 check, ad ognuno scelgo quello più vicino al cerchio unitario
for i=1:M
    [val,ind]=min(abs(abs(rad_svd)-1));
    w_svd(gg,i)=-angle(rad_svd(ind));
    rad_svd=[rad_svd(1:ind-1);rad_svd(ind+1:end)];
end

%disp('Soluzioni con inversione SVD e problema sovradeterminato: ');
w_svd(gg,:);
w_svd(gg,:)=sort(w_svd(gg,:));
%end
end

w_svd_tot=w_svd_tot+w_svd;
w_LS_tot=w_LS_tot+w_LS;
end

w_svd_tot=w_svd_tot./Mitr;
w_LS_tot=w_LS_tot./Mitr;


%Calcolo MSE

%MSE_QR=abs(w_svd_tot.^2-repmat(vere,size(w_svd_tot,1),1).^2);
MSE_LS=abs(w_LS_tot.^2-repmat(vere,size(w_LS_tot,1),1).^2);
MSE_SVD=abs(w_svd_tot.^2-repmat(vere,size(w_svd_tot,1),1).^2);




figure,...
    plot(SNR,MSE_LS), title('MSE della stima LS per ogni frequenza'),...
    xlabel('SNR [dB]'),ylabel('MSE'),...
    legend('puls norm 1','puls norm 2','puls norm 3'),...
    figure,plot(SNR,sum(MSE_LS'),'k'),title('MSE totale della stima LS di frequenza con HOYW'),...
    xlabel('SNR [dB]'),ylabel('MSE'),...
    
    
    subplot(2,1,2),...
    plot(SNR,MSE_SVD), title('MSE della stima SVD per ogni frequenza'),...
    xlabel('SNR [dB]'),ylabel('MSE'),...
    legend('puls norm 1','puls norm 2','puls norm 3'),...
%     subplot(3,1,3),...
%     plot(SNR,MSE_SVD), title('MSE della stima SVD per ogni frequenza'),...
%     xlabel('SNR [dB]'),ylabel('Pulsazione normalizzata [rad]'),...
%     


