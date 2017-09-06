%prova filtraggio lineare con matrice di convoluzione

clear all;
close all;
clc;

fc=1000;
t=linspace(0,1,fc);

y=sin(2*pi*100*t);

h=[1 -0.2];

n0=30;
n=5;

U=convmat(y,h,n0,n);

%genero le osservazioni

w=randn(n,1);

x=U*h'+w;

%calcolo DFT

W=dft(y',length(y));

Y=W*y';

fk=assefr(length(y),fc);

figure,plot(fk,abs(fftshift(Y)))
