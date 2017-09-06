%% Trasformata Z 

clear all;
close all;
clc;

%Partendo dalla definizione del filtro sequenza (N+1 campioni --> N zeri --> N sequenze elementari di 2 campioni)

a=[1 -0.5-j*0.2 1];

%trasf Z

%Voglio valutare la funzione di traferimento della sequenza, ovvero come
%questa cambia in fase e ampiezza una sinudoide complessa campionata a frequenza fc,
%in ingresso, per ogni frequenza da 0 a fN=fc/2

fc=2e3; %freq di campionamento delle sinusoidi complesse campionate usate per il calcolo 
%della FDT
fN=fc/2;

%Potrei a questo punto valutare la fdt per frequenze continue, z=exp(j*2*pi*(f/fc))


M=300; %numero di frequenze generalizzate nelle quali valuto la FDT
f=-(fc/2)+fc/M:fc/M:fc/2;
for k=1:length(f)
A=1;
for i=2:length(a)
    A=A+a(i)*exp(-j*2*pi*(f(k)/fc));
end
B(k)=A;
clear A;
end

figure,subplot(1,2,1),plot(f,abs(B)),title('Modulo della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Ampiezza'),axis([-fN fN 0 max(abs(B))+1]),...
    subplot(1,2,2),plot(f,angle(B)),title('Fase della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Fase [rad]'),axis([-fN fN -pi pi]),...

%% Trasformata Z: esempio da Piano Z --> sequenza || creazione filtri

%Ora al contrario: parto dalla definizione in frequenza di come voglio il
%filtro, cioè dal posizionamento degli zeri, ed arrivo all'espressione
%della sequenza corrispondente al filtro desiderato
clear all;
close all;
clc;

%Parto da un numero di zeri che è N, otterrò una sequenza che sarà
%l'antitrasformata Z del polinomio caratterizzato da questi zeri; ovvero la
%convoluzione delle sequenze elementari formate da 2 campioni: il primo
%pari a 1, ed il secondo pari a -z(k) ---> sequenza da (N+1) campioni


z=[0.9 0.5+0.5*1i 2*1i]; %zeri nel piano Z --> azzero la freq costante e fN/2

%Per avere una funzione di trasferimento (risposta in frequenza) a fase minima basta 
%cercare di contenere tutti gli zeri dentro al cerchio unitario 

%Uno o + zeri posti fuori dal cerchio unitario, se ne esiste anche qualcuno
%dentro, originano una risposta in frequenza a fase mista, nella quale si
%nota che, se lo zero a fase massima fosse uno solo, più questo viene
%spinto fuori dal cerchio unitario, lungo la sua frequenza di pertinenza,
%più il salto di fase corrispondente alla sua presenza shifta in frequenza
%verso frequenze + alte di quella di sua pertinenza!

%Traccio la funzione di trasferimento corrispondente ad un filtro con
%trasformata Z che ha questi zeri

fc=2e3; %freq di campionamento delle sinusoidi complesse campionate usate per il calcolo 
%della FDT
fN=fc/2;

%Potrei a questo punto valutare la fdt per frequenze continue, z=exp(j*2*pi*(f/fc))


M=300; %numero di frequenze generalizzate nelle quali valuto la FDT
f=-(fc/2)+fc/M:fc/M:fc/2;
for k=1:length(f)
A=1;
for i=1:length(z)
    A=A*(1-z(i)*exp(-1i*2*pi*(f(k)/fc)));
end
B(k)=A;
clear A;
end

figure,subplot(1,2,1),plot(f,10*log10(abs(B))),title('Modulo della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Ampiezza [dB]'),axis([-fN fN min(10*log10(abs(B)))-1 max(10*log10(abs(B)))+1]),...
    subplot(1,2,2),plot(f,angle(B)),title('Fase della FDT fino a fN'),...
    xlabel('f [Hz]'),ylabel('Fase [rad]'),axis([-fN fN -pi pi]),...

%Perfetto: ora voglio ottenere la sequenza

A=[1 -z(1)];
for i=2:length(z)
    A=conv(A,[1 -z(i)]);
end

n=0:length(A)-1;

figure,subplot(1,2,1),stem(n,real(A)),xlabel('indice [n]'),ylabel('Ampiezza'),title('Sequenza del filtro - parte reale'),...
    subplot(1,2,2),stem(n,imag(A)),xlabel('indice [n]'),ylabel('Ampiezza'),title('Sequenza del filtro - parte immaginaria'),
    


