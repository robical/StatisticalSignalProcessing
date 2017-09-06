%% Test di stima di potenza campionaria ddi un processo gaussiano incor

clear all;
close all;
clc;

N=1635;

sig=3; %deviazione standard!! non è la varianza, che in questo caso è 9, dato che il valore medio è 0

%for i=1:100 %100 realizzazioni del processo, plotto su un grafico l'errore
x=sig*randn(1,N); %1 realizzazione del processo

S_st=1/N*(sum(x.^2)); %stimatore non polarizzato ottimo

%err(i)=abs(S_st-sig^2);
%end

%media_er=mean(err);

disp(strcat('La potenza del processo vera è: ',num2str(sig^2)));
disp(strcat('La potenza stimata è: ',num2str(S_st)));

%figure,plot(1:100,err,'k'), xlabel('indice Realizzazione'),ylabel('Modulo errore stima'),...
%    title('Stima di potenza di processo gaussiano incorrelato'),hold on,
%plot(1:100,media_er,'r*')

%Utile:
% Generate values from a bivariate normal distribution with specified mean
%        vector and covariance matrix.
%           mu = [1 2];
%           Sigma = [1 .5; .5 2]; %Matrice di covarianza desiderata
%           R = chol(Sigma); %decompongo con Cholesky, che fattorizza in
%           una unica matrice radice, la matrice Sigma=R^2
%           z = repmat(mu,100,1) + randn(100,2)*R;

