%% Multislot estimation of frequency selective fast varying channel
% Alias:  Progettino Rampino

clear all;
close all;
clc;

%Caratteristiche del sistema di tx

Rs=100; %100 simboli/s
rho=0; %roll off - 0 è banda MIN
B=Rs*(1+rho); %banda occupata
Ts=1/B; %tempo di simbolo
Nfft=256;
dur=Ts;
%numero di campioni per ogni simbolo dispari - centro il campione centrale
%della forma d'onda nell'istante Ts
N=32; %campioni della forma d'onda
fc=N/Ts; %frequenza di campionamento della forma d'onda (per simulazione con Matlab
% in realtà c'è la conversione ADC)

fk=assefr(Nfft,fc);

%Modulazione BPSK
Rb=Rs;
stream=randi([0 1],1,Rs); %flusso di simbolo, con Rate di simbolo Rs
stream(stream==0)=-1; %modulazione BPSK

%Creo filtro formatore

Ht=ones(1,N)*sqrt((1/N)); %energia unitaria

%segnale analitico
H_tx=hilbert(Ht); %equivalente passa basso

%sovracampionamento dell'asse temporale per ospitare la forma d'onda
stream_tx=upsample(stream,N);
flusso_sim=conv(stream_tx,H_tx);
%flusso_sim=flusso_sim(1:N*Rs);
figure,stem(flusso_sim)

%creo filtro canale - equivalente passa basso

%per non avere selettività in frequenza, il max ritardo relativo tra i
%cammini in un canale multipath deve essere < del tempo di simbolo;  in
%questo caso non ho distorsione, ed il canale è rappresentabile attraverso
%un coefficiente complesso (come eq passa basso), dato che la forma d'onda
%ricevuta non cambia; tuttavia i ritardi sono trascurabili solo in valore
%assoluto, e non per il loro effetto sulla fase, che invece deve essere
%tenuto in conto specie se si lavoro in alta frequenza

%per un modello a 2 raggi, potrei pensare che il canale è tempo variante da
%simbolo a simbolo perchè il terminale si muove, e dunque il ritardo varia
% --> ho fatto un modello a 2 raggi completo che contiene anche tempo
% varianza del canale associata al movimento del terminale, potremmo usare
% quello per il momento (?)

d0= 2000; %distanza in metri iniziale MS-BS (a livello del terreno)
d=d0;
c=3*1e8; %velocità di propagazione in [m/s]
hBS=24; %altezza della BS [m] 
hMS=1.7; %altezza della MS [m] - uomo che parla e passeggia/corre
v=15/3.6; %velocità dell'uomo che cammina a 5 Km/h -> in [m/s] basta dividere per 3.6
dh=hBS-hMS; %differenza altezza antenne BS-MS

tau0=sqrt(dh^2+d0^2)/c; %ritardo iniziale di propagazione del raggio diretto

%ritardo di propagazione raggio riflesso
[l1,l2]=fermat(hBS,hMS,d); %funzione che trova i cammini minimi BS-MS
taur=(l1+l2)/c; %ritardo iniziale di propagazione raggio riflesso

dtau=taur-tau0; %prendo sempre il ritardo del raggio diretto come riferimento,
%uso solo il ritardo relativo

a=-0.1;
e=0.02;

for i=1:Rs
    % a = numero reale che esprime l'ampiezza in modulo associata al raggio
    % riflesso - cambia in teoria ad ogni Ts, perchè cambia l'angolo di
    % incidenza del raggio riflesso sulla superficie (terreno)
    % -al momento lo tengo fisso- 
    a=a-e; %riflessione totale, viviamo su uno specchio
    
    %modello del canale passa basso
    h(:,i)=1+a*exp(-1i*2*pi*fk*dtau); %modello a 2 raggi

    %coe(i)=sum(h); %coefficiente complesso rappresentante il canale, tempo invariante ALL'INTERNO
    %DEL TEMPO DI SIMBOLO --> non selettivo in frequenza
    %dunque se si volesse stimare il canale su una sequenza di training
    %composta da N simboli, troverei N coefficienti complessi rappresentativi
    %del canale, uno per tempo di simbolo

    %usando il modello MS-BS che ho fatto, la dipendenza è:
    deltad=v*Ts; %distanza che percorre omino che cammina per strada 
    
    %per raggio diretto
    dthe=deltad*(1/(1+(d/dh)^2));
    dl=dh/(cos(dthe)); %variazione relativa del percorso diretto di propagazione
    tau0=tau0+(dl/c); %nuovo ritardo di propagazione per cammino diretto
    
    d=d+deltad; %nuova posizione della MS
    
    %per raggio riflesso
    [l1,l2]=fermat(hBS,hMS,d); %funzione che trova i cammini minimi BS-MS, per la nuova posizione della MS
    taur=(l1+l2)/c; %ritardo iniziale di propagazione raggio riflesso
    
    %ritardo relativo
    dtau=taur-tau0;
end

    %filtro di canale tempo variante
    %h_f=[ones(size(h));h];
    
    
    %Qua fa fatta la convoluzione in tempo o moltiplicazione in frequenza
    %per il canale
    
    
    temp=zeros(size(flusso_sim));
    N_conv=N+size(h,1)-1; %dim conv aperiodica
    Ht_pad=[H_tx zeros(1,N_conv-length(H_tx))];
    T2=fftshift(fft(Ht_pad,Nfft));
    for i=1:Rs
        %T1(:,i)=fft(h_f_pad,Nfft);
        T3=h(:,i)'.*T2;
        temp((i-1)*N+1:i*N+1)=temp((i-1)*N+1:i*N+1)+(1/Nfft)*ifft(fftshift(T3),N+1); %convoluzione del simbolo con il filtro al tempo di simbolo 
        %corrispondente
    end
        
    %Qui ci va il filtro inverso, stimato o ottenuto in qualche modo
    new_temp=zeros(size(flusso_sim));
    for i=1:Rs
        slice=temp((i-1)*N+1:i*N+1); %è già della lunghezza della conv aperiodica
        Tslice=fftshift(fft(slice,Nfft));
        Tdec=(1./h(:,i))'.*Tslice; %deconvoluzione
        new_temp((i-1)*N+1:i*N+1)=new_temp((i-1)*N+1:i*N+1)+(1/Nfft)*ifft(fftshift(Tdec),N+1);
    end
        
    
    %filtro adattato in ricezione
    H_ad=fliplr(conj(H_tx));
    
    %convoluzione con filtro adattato (sul simbolo), e campionamento nel
    %MAX
    
    sig_rx=zeros(size(flusso_sim));
    for i=1:Rs
        slice=new_temp((i-1)*N+1:i*N+1); %è già della lunghezza della conv aperiodica
        Tslice=fftshift(fft(slice,Nfft));
        Tdec=(1./h(:,i))'.*Tslice; %deconvoluzione
        sig_rx((i-1)*N+1:i*N+1)=sig_rx((i-1)*N+1:i*N+1)+(1/Nfft)*ifft(fftshift(Tdec),N+1);
    end
    %la parte immaginaria è dovuta a errori numerici, la elimino a mano
%     s_i=imag(sig_rx);
%     sig_rx=sig_rx-1i*s_i; %fatto

    %figure,plot(1:length(sig_rx),sig_rx) %segnale ricevuto non campionato

    %campionamento a Rs, ricordare che essendo i filtri causali, c'è un ritardo
    %totale di tau0=Ts, dal quale si comincia a campionare a k*Ts, Ts/2 dovuto
    %al fitro di TX e Ts/2 dovuto al filtro adattato
    sig_rx_camp=downsample(sig_rx(N:end-N),N);
    
    
    figure,subplot(3,2,1),stem((1:Rs)*Ts,stream),xlabel('Tempo [n*Ts] [s]'),ylabel('Ampiezza'),title('Flusso di simboli inviati'),hold on,...
        subplot(3,2,3),stem((1:Rs)*Ts,real(sig_rx_camp),'b'),xlabel('Tempo [n*Ts] [s]'),ylabel('Ampiezza'),title('Flusso di simboli ricevuti - reale'),hold on,...
        subplot(3,2,4),stem((1:Rs)*Ts,imag(sig_rx_camp),'r'),xlabel('Tempo [n*Ts] [s]'),ylabel('Ampiezza'),title('Flusso di simboli ricevuti - imag'),hold on,...
        subplot(3,2,5),stem((1:Rs)*Ts,abs(sig_rx_camp),'g'),xlabel('Tempo [n*Ts] [s]'),ylabel('Ampiezza'),title('Flusso di simboli ricevuti - modulo'),hold on,...
        subplot(3,2,6),stem((1:Rs)*Ts,angle(sig_rx_camp),'m'),xlabel('Tempo [n*Ts] [s]'),ylabel('Fase'),title('Flusso di simboli ricevuti - fase'),