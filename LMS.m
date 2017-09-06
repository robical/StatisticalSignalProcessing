%INPUT
%
%
% U= matrice di convoluzione per identificazione
% h_sti= stima iniziale del filtro
% y= osservazioni con rumore
% mu= passo di aggiornamento
% h=filtro vero
%
%
%OUTPUT:
%
%
% fstim=filtro stimato
% MSE=errore quadratico medio della stima ad ogni iterazione


function [fstim,MSE]=LMS(U,h_sti,y,mu,h)


P=length(h_sti);

i=1;
while(i==1 || i==2 || J(i-1)<J(i-2) || (abs(J(i-1)/J(i-2))>(1-1e-2) && abs(J(i-1)-J(i-2))>1e-3))
    %stima dell'out
    
    the(i)=U(i,:)*h_sti(:,i); %scalare, stima dell'OUT al passo i
    
    %calcolo del funzionale di costo
    
    %stima campionaria della crosscorrelaz istantanea
    p=(y(i)'*U(i,:)); %riga
    
    %stima campionaria della matrice di covarianza istantanea
    R=zeros(P);
    R=(U(i,:)'*U(i,:)); 
    
    J(i)=(y(i)'*y(i))-2*p*h_sti(:,i)+h_sti(:,i)'*R*h_sti(:,i);
    
    %aggiornamento filtro
    
    err(i)=y(i)-U(i,:)*h_sti(:,i); %errore istantaneo di stima
    h_sti(:,i+1)=h_sti(:,i)+mu*(err(i)*U(i,:)');
    i=i+1;
end

fstim=h_sti(:,end); %stima del filtro a convergenza

%errore quadratico medio di stima ad ogni iterazione
err=h_sti-repmat(h',1,size(h_sti,2));

for i=1:size(err,2)
    MSE(i)=err(:,i)'*err(:,i); %calcolo MSE tra filtro vero e filtro stimato alle varie iterazioni
end