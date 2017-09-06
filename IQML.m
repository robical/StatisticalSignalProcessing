%% IQML-QC (Iterative quadratic ML con vincolo quadratico)
%Stima di frequenza ML con metodo iterativo, 
%le ottimizzazioni ad ogni iterazione avvengono in 1 passo attraverso
%ottimizzazione vincolata di una forma quadratica, grazie alla
%parametrizzazione del problema come la ricerca dei parametri di un filtro
%che abbia poli posizionati esattamente (o il più vicino possibile) ai poli
%del processo, che contenendo M sinusoidi (pure per il momento), ha M poli
%sul cerchio unitario

%Funziona solo se non c'è rumore!! ORA FUNZIONA ANCHE CON RUMORE!!!!

clear all;
close all;
clc;

%Creo le osservazioni

M=3; %3 sinusoidi pure

N=100; %lunghezza osservazione
%sign=1;

sign=2; %Potenza totale del rumore gaussiano circolare complesso
Amp=1;
SNR=20*log10(Amp/sign);

%Fasi di ogni sinusoide
phi1=0.1;
phi2=0.2;
phi3=0.3;


Mitr=100; %numero massimo iterazioni

disp(strcat('L''SNR delle osservazioni è: ',num2str(SNR),'dB'));
%Nfft=2^12;

%chiaramente devo definire pulsazioni normalizzate comprese tra -pi e pi
w1=pi/2;
w2=pi/4;
w3=pi*(2/3);

%Devo modificare un po' il problema per tenere conto di fasi e ampiezza
%casuali

%matrice A(w)

A=[exp(1i*w1*[0:N-1]') exp(1i*w2*[0:N-1]') exp(1i*w3*[0:N-1]')];

%vettore di fasi e ampiezze casuali, va nel ciclo,lo chiamo s
% fasi --> distibuite in modo uniforme tra 0 e 2*pi
%ampiezze --> distibuite Rayleigh


%Inizializzo la matrice di convoluzione fatta con l'ipotetico filtro

itr=1; %counter delle iterazioni

b(1,:)=ones(1,M+1); %inizializzo la soluzione

%Creazione del modello delle osservazioni



%Devo scegliere come voglio crearlo, si potrebbe fare che per ogni istante
%temporale N, il rumore cambia, le fasi cambiano, ma le ampiezze restano
%costanti nella finestra che poi utilizzo per la stima

%Caso con : ampiezze fisse (det) e fasi e rumore casuale in ogni istante


s=repmat(Amp,M,1).*exp(1i*[phi1;phi2;phi3]); %vettore delle ampiezze det e fasi det - raylrnd()
    
    
 %va eliminato l'errore derivante da sfasamenti fissi
%  x=x_fas(2:end)./x_fas(1:end-1); %Raddoppia la potenza di rumore farlo
%  così!!
%  x=[x_fas(1) x];


x_sum=zeros(1,N); %realizzazioni mediate, le uso per stimare poi LS ampiezze e fasi

for itr=2:Mitr
    
%nuova realizzazione di rumore ad ogni iterazione
x=(A*s)'+(sqrt(sign/2))*(randn(1,N)+1i*randn(1,N));
x_sum=x_sum+x;
%Creo matrice di convoluzione
%in blocchi da M=3
k=1;
for i=1:1:N-M
    X(k,:)=x(i:i+M);
    k=k+1;
end

k=1;
for i=1:1:N-M
    B(k,:,itr)=[zeros(1,i-1) fliplr(b(itr-1,:)) zeros(1,N-M-(i-1))];
    k=k+1;
end



%matrice della forma quadratica da ottimizzare
%Uso QR
temp=B(:,:,itr)*B(:,:,itr)';
[Q,R]=qr(temp);

C(:,:,itr)=X'*(R\Q')*X;

%Vincolo quadratico

[V,D]=eig(C(:,:,itr)); %trovo autovalori e autovettori, scegliere il più piccolo
%equivale a imporre il vincolo di norma 2 del vettore b ==1

[y,ind]=min(diag(D));

b(itr,:)=V(:,ind)'; %stima dei coefficienti (M+1)

%devo trovare come sono disposti i coefficienti
poli(:,itr)=roots(fliplr(b(itr,:))); %proviamo così

%le pulsazioni sono le fasi dei poli

puls_norm(:,itr)=abs(angle(poli(:,itr)));
end

x_sum=x_sum/Mitr; %realizzazioni mediate

%Andamento di polarizzazione ed MSE dello stimatore
w_vere=[w2;w1;w3];

%Rimuovo qui la dipendenza dalle fasi

% nuove=2*puls_norm(:,2:1:end)-puls_norm(:,1:1:end-1);
% puls_norm=[puls_norm(:,1) nuove];

for itr=1:Mitr
    tem=zeros(M,1);
    for k=1:itr
        tem=tem+puls_norm(:,itr);
    end
        tem_pol(:,itr)=tem./itr - w_vere; %polarizzazione della stima di ciascuna frequenza separata
        pol(itr)=sum(tem_pol(:,itr)); %polarizzazione totale al variare delle iterazioni
end

figure,plot(1:itr,tem_pol),title(strcat('Andamento polarizzazione stimatore all''aumentare delle iterazioni con SNR ',num2str(SNR),'dB')),...
    xlabel('# iterazioni'),ylabel('Polarizzazione'),...
    legend('Puls norm 2','Puls norm 1','Puls norm 3')

%Andamento MSE all'aumentare delle iterazioni 

for itr=1:Mitr
    tem=zeros(M,1);
    for k=1:itr
        tem=tem+(abs(puls_norm(:,itr).^2-w_vere.^2)); %varianza
    end
        MSE_single(:,itr)=tem./itr + abs(tem_pol(:,itr)).^2; %polarizzazione al variare delle iterazioni
        MSE_tot(itr)=sum(tem)/itr + abs(pol(itr))^2;
end

figure,plot(1:itr,MSE_single),title(strcat('Andamento MSE stimatore all''aumentare delle iterazioni con SNR ',num2str(SNR),'dB')),...
    xlabel('# iterazioni'),ylabel('MSE'),...
    legend('Puls norm 2','Puls norm 1','Puls norm 3')

%Recupero i valori di fase e ampiezza

A=[exp(1i*puls_norm(2,end)*[0:N-1]') exp(1i*puls_norm(1,end)*[0:N-1]') exp(1i*puls_norm(3,end)*[0:N-1]')];

%Stima LS del vettore con ampiezze e fasi

s_stim=(A'*A)\(A'*x_sum');

amp=abs(s_stim);
fasi=angle(s_stim);

disp('Le ampiezze stimate LS: ');
amp
disp('Le fasi stimate LS: ');
fasi

