%% prova

clear all;
close all;
clc;


N=4;

W=exp(1i*((2*pi)/N));

for k=0:(N-1)
    for n=0:(N-1)
        Wd(k+1,n+1)=W^(-n*k); %trasf diretta DFT
    end
end

Wi=(1/N)*Wd'; %trasf DFT inversa


%La differenza sta nella distribuzione della potenza tra parte reale ed
%immaginaria , tra i 2 processi reale e complesso
% La DFT diagonalizza anche il processo filtrato, e mostra dunque la stima
% della potenza alle frequenze discrete rappresentate dagli autovettori
% della matrice di covarianza

%Creazione filtro

h=[1 1]'; %filtro da 2 campioni (colonna)
M=length(h); %lunghezza del filtro
h=[h' zeros(1,N-M)]'; %filtro con zero padding (colonna)
hl=[h' h']; %vettore che contiene 2 periodi di h (riga) 

%creazione matrice di convoluzione circolare
for i=0:(N-1)
    H(i+1,:)=hl(end-N-i+1:end-i);
end


%Inizializzo matrice di covarianza per la stima, a zero

RX=zeros(N,N);
Rx=zeros(N,N);
Ry=zeros(N,N);
RY=zeros(N,N);
RXc=zeros(N,N);
Rxc=zeros(N,N);

L=200000; %numero di realizzazioni del processo casuale usate per la stima della matrice di covarianza

for l=1:L
%segnale bianco reale

x=randn(N,1);

Rx=(x*x')+Rx;

%convoluzione circolare della sequenza reale bianca con il filtro

y=H*x;

Ry=(y*y')+Ry;

% %segnale bianco complesso
% xc=(1/sqrt(2))*randn(N,1)+1i*(1/sqrt(2))*randn(N,1); %parte reale e immaginaria si sommano in potenza!!
% 
% Rxc=(xc*xc')+Rxc;


%prova che le componenti della sua DFT sono indipendenti (segnale reale)

X=Wd*x;

%Stima della matrice di covarianza della DFT

RX=(X*X')+RX;

%DFT del segnale filtrato

Y=Wd*y;

RY=(Y*Y')+RY;

% %prova che le componenti della sua DFT sono indipendenti (segnale complesso)
% 
% Xc=Wd*xc;
% 
% %Stima della matrice di covarianza della DFT
% 
% RXc=(Xc*Xc')+RXc;

end

Rx=(1/L)*Rx; %stima matrice di covarianza del processo reale
RX=(1/L)*RX; %stima covarianza della DFT del processo reale

Ry=(1/L)*Ry; %stima matrice di covarianza del processo reale bianco filtrato (colorato)
RY=(1/L)*RY; %stima della matrice di covarianza della DFT del processo filtrato

Rxc=(1/L)*Rxc; %stima covarianza del processo complesso
RXc=(1/L)*RXc; %stima della covarianza della DFT del processo complesso

% disp('Stima della matrice di covarianza della DFT del processo (complesso)'),
% RXc
% disp('Stima della matrice di covarianza del processo (complesso - var re 1/sqrt(2) e im 1/sqrt(2))'),
% Rxc

disp('Stima della matrice di covarianza della DFT del processo (reale)'),
RX
disp('Stima della matrice di covarianza del processo (reale - var 1)'),
Rx

disp('Stima della matrice di covarianza della DFT del processo (reale) FILTRATO'),
RY
disp('Stima della matrice di covarianza del processo (reale - var 1) FILTRATO'),
Ry






