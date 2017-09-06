%Esempio di stimatore BLUE per modelli lineari

clear all;
close all;
clc;


%filtraggio lineare + rumore

f0=200;
fc=500;
t=linspace(0,100*1/fc-1/fc,fc);
y=sin(2*pi*f0*t);

n0=50;
n=100;

h=[1;-0.2];

U=convmat(y,h,n0,n);

%creo osservazioni
w=4*randn(size(U,1),1);

x=U*h+w;

figure,plot(x),title('Osservazioni rumorose')

Cw=4*eye(size(U,1)); %stima della matrice di covarianza del rumore

%stimatore BLUE

the=(U'*Cw^-1*U)^-1*U'*Cw^-1*x;

cov_the=the*the';

disp(strcat('La stima della matrice di covarianza di rumore che è stata usata è: '))

imagesc(Cw)

disp(strcat('La stima dei parametri del filtro è: the1 ',num2str(the(1))))
disp(strcat('La stima dei parametri del filtro è: the2 ',num2str(the(2))))

disp(strcat('La matrice di covarianza della stima è: '))

cov_the

%Ok dinchè conosco completamente! la matrice di covarianza del rumore,
%altrimenti sbaglia, anche parecchio. Per esempio un grande errore viene
%commesso considerando tutte le osservazioni (i campioni della sequenza rumorosi)
% incorrelati tra loro --> matrice di covarianza di rumore diagonale e
% dunque anche matrice di covarianza del segnale diagonale, anche
% ammettendo di conoscere la varianza delle singole variabili casuali
% componenti il vettore aleatorio che è l'osservazione di rumore (processo)