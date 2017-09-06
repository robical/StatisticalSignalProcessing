%% EKF per localizzazione e tracking (Esempio 3 BS e 1 MS)

clear all;
close all;
clc;

%Parametri
T=1; %campionamento della posizione della MS, 1 sec

%Posizioni delle 3 BS - fisse
rx1=10;
ry1=10;
rx2=120;
ry2=40;
rx3=56;
ry3=22;

gamma=100; %parametro legato allo spettro della forma d'onda trasmessa
c=3*1e8; %velocità di propagazione
Pt=100; %potenza costante trasmessa dalla MS [mW]
R0=10; %distanza di riferimento [m]
SNR0=100; %SNR alla distanza di riferimento

%Parametri per equazione dinamica dello stato
siga=1;
Caa=[0 0 0 0;0 0 0 0;0 0 siga 0;0 0 0 siga];
A=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];

Rbs=[rx1 ry1;rx2 ry2;rx3 ry3]; %posizioni delle BS
N=size(Rbs,1);

init=[5 5 2 2]; %inizializzazione stato MS

%vettore di stato the(n)=[rx(n) ry(n) vx(n) vy(n)]; contiene posizione 
%nel piano e velocità della MS

the_po(:,1)=init'; %inizializzazione della stima MMSE a 
the(:,1)=the_po(:,1);
%posteriori dello stato
P=size(the,1); %cardinalità spazio dei parametri


%inizializzazione della covarianza a posteriori
P_po(:,:,1)=eye(P);

for i=2:10

%Equazione dinamica dello stato

a=[0 0 sqrt(siga)*randn sqrt(siga)*randn]'; %vettore delle accelerazioni

the(:,i)=A*the(:,i-1)+a; %aggiornamento dello stato via equazione dinamica

%PREDIZIONE
the_pr(:,i)=A*the_po(:,i-1); %aggiornamento del valore medio a priori dello stato,
%cioè la predizione all'istante n date (n-1) osservazioni
P_pr(:,:,i)=A*P_po(:,:,i-1)*A'+Caa; %aggiornamento della covarianza a priori 
%cioè della predizione


%Legame stato-osservazioni
%Legame non lineare!
%le potenze del rumore dipendono dalla distanza BS-MS
sigmab(:,i)=(((c^2)*gamma)/SNR0)*((sqrt((repmat(the_pr(1,i),N,1)-Rbs(:,1)).^2+(repmat(the_pr(2,i),N,1)-Rbs(:,2)).^2))./(R0^2)); %covarianze di rumore
Cbb(:,:,i)=sigmab(:,i)*sigmab(:,i)'; %matrice di covarianza del disturbo, tempo variante
%tempo varianti
b(:,i)=sqrt(sigmab(:,i)).*randn(N,1);
R(:,i)=sqrt((the(1,i)-Rbs(:,1)).^2+(the(2,i)-Rbs(:,2)).^2)+b(:,i); %osservazioni

%AGGIORNAMENTO
%Linearizzazione della funzione non lineare
B(:,:,i)=[the_pr(1,i)/(sqrt((the_pr(1,i)-Rbs(1,1)).^2+(the_pr(2,i)-Rbs(1,2)).^2)) the_pr(2,i)/(sqrt((the_pr(1,i)-Rbs(1,1)).^2+(the_pr(2,i)-Rbs(1,2)).^2)) 0 0;...
    the_pr(1,i)/(sqrt((the_pr(1,i)-Rbs(2,1)).^2+(the_pr(2,i)-Rbs(2,2)).^2)) the_pr(2,i)/(sqrt((the_pr(1,i)-Rbs(2,1)).^2+(the_pr(2,i)-Rbs(2,2)).^2)) 0 0;...
    the_pr(1,i)/(sqrt((the_pr(1,i)-Rbs(3,1)).^2+(the_pr(2,i)-Rbs(3,2)).^2)) the_pr(2,i)/(sqrt((the_pr(1,i)-Rbs(3,1)).^2+(the_pr(2,i)-Rbs(3,2)).^2)) 0 0];

%Variabili d'appoggio
Ppr=P_pr(:,:,i);
Bb=B(:,:,i);
Cb=Cbb(:,:,i);
G(:,:,i)=(Ppr*Bb')/(Bb*Ppr*Bb'+Cb); %guadagno di kalman  tempo variante
the_po(:,i)=the_pr(:,i)+G(:,:,i)*(R(:,i)-sqrt((the_pr(1,i)-Rbs(:,1)).^2+(the_pr(2,i)-Rbs(:,2)).^2)); %aggiornamento stima MMSE (a posteriori)
%dello stato

%Equazione di Riccati
P_po(:,:,i)=P_pr(:,:,i)-G(:,:,i)*B(:,:,i)*P_pr(:,:,i); %aggiornamento covarianza della stima

end

%FAtto è fatto, ora bisogna capire perchè si impalla...