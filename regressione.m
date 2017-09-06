% Regressione lineare VS regressione polinomiale
clear all;
close all;
clc;


fc=1000;

t=0:1/fc:1-1/fc;
w=randn(size(t));



y=5+4*t+w;

wrid=w(1:2:end);
trid=t(1:2:end); %set ridotto di misure
yrid=y(1:2:end);

Hlin=[ones(size(trid))' trid'];
Hpol=[ones(size(trid))' trid' (trid.^2)'];


figure,plot(trid,yrid,'r'),title('Stima LS con modello lineare (nero) e polinomiale-grado2 (verde)')

%soluzione ai minimi quadrati con modello lineare

theta_s=(Hlin'*Hlin)^-1*Hlin'*yrid';

y_ric=theta_s(1)+theta_s(2)*trid;

hold on, plot(trid,y_ric,'k'),

%soluzione ai minimi quadrati con modello polinomiale (2 grado)

theta_s=(Hpol'*Hpol)^-1*Hpol'*yrid';

y_ric=theta_s(1)+theta_s(2)*trid+theta_s(3)*trid.^2;

hold on, plot(trid,y_ric,'g'),