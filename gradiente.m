%% Algoritmo LMS applicato all'identificazione di un filtro a 3 prese

clear all;
close all;
clc;


% creo filtro

e=1e-2;
A=[(1-e)*exp(1i*(pi/6)) (1-e)*exp(-1i*(pi/6))];

[h]=zero(A); %crea filtro FIR con la sequenza di zeri che gli ho dato, 
%causale e stabile (è un passa alto un po' scarso)

fc=2*1e3;
fN=fc/2;

[Spettro,f]=spet_plo(A,[0],fc);


%Ora serve un ingresso bello lungo

Tc=(1/fc);

dur=4; %durata = 3 secondi

t=0:Tc:dur-Tc;

%INGRESSO: processo gaussiano (K=3) bianco stazionario di potenza sigx

sigx=5;

x=sqrt(sigx)*randn(1,length(t));

%Ora da questo creo la matrice di convoluzione
%P=3
%N=6000+3-1=6002
%Verrà una matrice di convoluzione 6002x3 !!

P=length(h); %lunghezza del filtro
N=length(x);

x=[ones(1,P-1)*0.5 x]; %meglio inizializzarlo con valori diversi da zero, 
% % altrimenti la prima stima della matrice di covarianza non permette di proseguire 

for i=1:N
    U(i,:)=x(i:1:i+P-1);
end

%Ora posso fare la convoluzione!

sigw=1;

w=sqrt(sigw)*randn(N,1); %vettore rumore gaussiano bianco potenza sigw

y=U*h'+w; %ecco le osservazioni rumorose

%Ora suppongo di essere a conoscenza della sequenza in ingresso, (sequenza di training)
% e provo a stimare con LMS i coefficienti del filtro (incognita), usando
% l'OUT più rumore, cioè y

%uso scelta conservativa per garantire la convergenza passo=1/traccia;

%scelta iniziale del filtro

h_sti=ones(length(h),1); %colonna

%da qui parte l'algoritmo LMS
%comincio senza graadiente stocastico

i=1;
tic;
while(i==1 || i==2 || J(i-1)<J(i-2) || (abs(J(i-1)/J(i-2))>(1-1e-6) && abs(J(i-1)-J(i-2))>1e-6))
    %stima dell'out
    
    the(i)=U(i,:)*h_sti(:,i); %scalare, stima dell'OUT al passo i
    
    %calcolo del funzionale di costo
    
    %stima campionaria della crosscorrelaz al passo i
    p=(1/i)*(y(1:i)'*U(1:i,:)); %riga
    
    %stima campionaria della matrice di covarianza
    R=zeros(P);
    R=(1/i)*(U(1:i,:)'*U(1:i,:)); %stima della mat covarianza con le prime i osservazioni
    
    J(i)=(1/i)*(y(1:i)'*y(1:i))-2*p*h_sti(:,i)+h_sti(:,i)'*R*h_sti(:,i);
    
    %aggiornamento filtro
    mu=(1/trace(R)); %scelta conservativa, non stimo gli autovalori di R --> convergenza più lenta ma garantita
    h_sti(:,i+1)=h_sti(:,i)+mu*(p'-R*h_sti(:,i));
    i=i+1;
end
toc;

%errore quadratico medio di stima ad ogni iterazione
err=h_sti-repmat(h',1,size(h_sti,2));

for i=1:size(err,2)
    MSE(i)=err(:,i)'*err(:,i); %calcolo MSE tra filtro vero e filtro stimato alle varie iterazioni
end


%Funziona!!! Ora proviamo con il gradiente stocastico

%Approssimo le correlazioni con il loro valore puntuale
mu=1/(3*trace(R));
h0=ones(length(h),1);
tic;
[fil_LMS,MSElms]=LMS(U,h0,y,mu,h);
toc;

%LMS normalizzato
h0=ones(length(h),1);
mu=0.04;
tic;
[fil_nLMS,MSEnlms]=nLMS(U,h0,y,mu,h);
toc;

%RLS
tic;
[fil_RLS,MSErls]=RLS(U,y,h);
toc;    

%Plotting dei risultati

figure,subplot(2,4,1),stem(h_sti(:,end)),xlabel('Prese del filtro'),ylabel('Ampiezza'),title('Filtro stimato col gradiente'),...
    subplot(2,4,2),stem(fil_LMS),xlabel('Prese del filtro'),ylabel('Ampiezza'),title('Filtro stimato con LMS'),...
    subplot(2,4,3),stem(fil_nLMS),xlabel('Prese del filtro'),ylabel('Ampiezza'),title('Filtro stimato con LMS normalizzato'),...
    subplot(2,4,4),stem(fil_RLS),xlabel('Prese del filtro'),ylabel('Ampiezza'),title('Filtro stimato con RLS'),...
    
    subplot(2,4,5),stem(MSE),xlabel('Iterazioni'),ylabel('Errore di stima'),title('MSE della stima con gradiente'),...
    subplot(2,4,6),stem(MSElms),xlabel('Iterazioni'),ylabel('Errore di stima'),title('MSE della stima con LMS'),...
    subplot(2,4,7),stem(MSEnlms),xlabel('Iterazioni'),ylabel('Errore di stima'),title('MSE della stima con LMS normalizzato'),...
    subplot(2,4,8),stem(MSErls),xlabel('Iterazioni'),ylabel('Errore di stima'),title('MSE della stima con RLS'),...
    
%Per misurare bene le prestazioni, potrebbe essere bello:
%- costruire un filtro di Wiener con il filtro stimato, e vedere qual'è
%quello che ottiene meglio il segnale in ingresso (sebbene a questo punto cada un po' il senso 
%di stimare i coefficienti del filtro, sarebbe bastato stimare l'ingresso)