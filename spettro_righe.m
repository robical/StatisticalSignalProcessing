%Analisi spettro a righe (DFT)

%Stima di fase e frequenza e ampiezza di una sinusoide immersa in rumore gaussiano
%bianco (dunque con vc estratte incorrelate), usando + realizzazioni del
%processo, cioè più acquisizioni

clear all;
close all;
clc;

%osservazione: somma sinusoide+rumore

N=100;  %numero di campioni componenti l'osservazione

fc=100;
fn=fc/2;
%t=linspace(0,(N-1)/fc,fc);


%sele=randi([1,fc],N,1);

%trid=t(sele);

%Modello

A=15; %ampiezza
phi=0.02; %fase
f0=10; %frequenza
w0=2*pi*f0; %cosinusoide a 20Hz


rum=1:0.5:12; %diversi valori della potenza di rumore

kk=1;
%for g=1:length(rum)
g=100:50:800;
for kk=1:length(g)    
    sigman=rum(4); %potenza di rumore
    SNR=20*log10(A/sigman); %SNR per osservazioni con rumore di potenza rum(g)
    t=linspace(0,(g(kk)-1)/fc,g(kk));



    l=10; %numero di realizzazioni (per abbassare l'effetto del rumore)
    tl=zeros(1,g(kk));
    for i=1:l
        n=sigman*randn(1,g(kk));
        x=A*cos(w0*t+phi)+n;
        tl=tl+x;
    end

    x=tl/l; %ora x è il vettore somma di tutte le realizzazioni della sinusoide + rumore
    x=x';
    
    [amp_stim(kk) freq_stim(kk) fas_stim(kk)]=stimAFF(x,t,fc); %Stima ML amp freq e fase  
%end

end


param=[amp_stim;freq_stim;fas_stim]; %vettore che contiene tutti i parametri stimati,
%per ogni valore della potenza di rumore scelto

figure,subplot(1,3,1),...
    plot(g,amp_stim,'r'),title('Stima ampiezza al variare di N'),...
    subplot(1,3,2),...
    plot(g,fas_stim,'b'),title('Stima fase al variare di N'),...
    subplot(1,3,3),...
    plot(g,freq_stim,'k'),title('Stima frequenza al variare di N'),...
    
