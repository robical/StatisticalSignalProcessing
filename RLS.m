%Algoritmo RLS per la stima dei coefficienti del filtro




function [fil,MSE]=RLS(U,y,h)


P=length(h);

delta=0.3;

R_s= eye(P)*delta; %inizializzazione della matrice di covarianza

f_s=zeros(P,1); %inizializzazione della stima

Ri=(R_s)^-1;

i=1;
while(i<=P+20*P)
    
    Ri=Ri - (Ri*U(i,:)'*U(i,:)*Ri)/(1+ U(i,:)*Ri*U(i,:)'); %aggiornamento dell'inversa
    
    g=Ri*U(i,:)'; %guadagno al passo i
    
    f_s(:,i+1)= f_s(:,i) + g* (y(i)-U(i,:)*f_s(:,i)); %aggiornamento filtro
    
    i=i+1;
end

fil=f_s(:,end);

%errore quadratico medio di stima ad ogni iterazione
err= f_s - repmat(h',1,size(f_s,2));

for i=1:size(err,2)
    MSE(i)=err(:,i)'*err(:,i); %calcolo MSE tra filtro vero e filtro stimato alle varie iterazioni
end