%% Algoritmo EM per la stima di ritardo/ampiezza (f. d'onda note)

clear all;
close all;
clc;

fc=100; %campionamento a 100Hz

t=0:1/fc:4; %3 secondi di osservazione

%fattore di forma della Ricker

f=3; %frequenza dominante della Ricker, 15hz

%Parametri

taumin=30/fc;

A1=6;
tau1=100/fc; %ritardo in secondi, conviene passare dai campioni e poi dividere per fc
A2=6;
tau2=130/fc; %ritardo in secondi, conviene passare dai campioni e poi dividere per fc
A3=6;
tau3=180/fc; %ritardo in secondi, conviene passare dai campioni e poi dividere per fc

Nstep=50; %numero di passi EM
sign=1; %potenza di rumore
costmin=10000; %inizializzazione

beta(1,:)=[1/3 1/3 1/3]; %pesa l'errore di predione (aggiornamento)

%stima iniziale dei parametri per l'avvio di EM

the(:,:)=[1 100/fc; 2 120/fc; 3 280/fc]; %contiene i parametri: Amp rit amp rit amp rit

rthe=[A1 tau1;A2 tau2;A3 tau3];

%forme d'onda

[w1,sig]=rickerw(t,f,A1,tau1);
[w2,sig]=rickerw(t,f,A2,tau2);
[w3,sig]=rickerw(t,f,A3,tau3);

%segnale somma di 3 Ricker con diversi ritardi e ampiezze


n=sqrt(sign)*randn(1,length(t));
Cn=sign*eye(length(t)); %covarianza di rumore

x = w1+w2+w3+n; %osservazioni incomplete, rumorose

%Ora assumo che le variabili complete siano i 3 segnali puliti con i loro
%parametri che dovrò stimare con EM, e su ognuno di questi segnali ci sarà
%un rumore con potenza 1/3, in modo che la somma in potenza equipari quella
%presente nelle osservazioni per ciò che attiene il rumore, dunque bk=1/3
%per ogni k

i=1;

aviobj = avifile('Cost_function1.avi');
aviobj2 = avifile('EM_evolve.avi');
 
fig=figure(1);
fig2=figure(2);

%griglia per verosimiglianza

A=0:0.5:10; %valori testati in ampiezza 
tau=taumin:5/fc:300/fc; %valori testati in ritardo
cost=zeros(length(A),length(tau),size(the,1));

for p=1:Nstep %solo 15 colpi per il momento
    
    % E-step
    
    %Come prima cosa creo segnale somma, cioè un modello delle
    %osservazioni, attraverso la prima stima dei parametri, che è
    %totalmente errata

        somma=zeros(3,length(t));
        somma(1,:)=rickerw(t,f,the(2,1),the(2,2))+rickerw(t,f,the(3,1),the(3,2));
        somma(2,:)=rickerw(t,f,the(1,1),the(1,2))+rickerw(t,f,the(3,1),the(3,2));
        somma(3,:)=rickerw(t,f,the(2,1),the(2,2))+rickerw(t,f,the(1,1),the(1,2));
%     somma=zeros(1,length(t));
%     for k=1:3
%         somma=somma+rickerw(t,f,the(k,1),the(k,2));
%     end
%     
    
    %poi stimo ogni segnale separatamente come somma di 2 componenti:
    % la prima: il modello del segnale singolo k, creato con i parametri
    % attuali
    % la seconda: l'errore di stima dell'osservabile, sempre con tutti i
    % parametri inseriti nel modello e con le osservazioni di cui dispongo,
    % pesando questo errore con la parte di rumore che attiene al segnale k
    
    for k=1:3
        y(:,k)=rickerw(t,f,the(k,1),the(k,2))+beta(k)*(x-somma(k,:)); %stima delle complete sulla base dei nuovi parametri
        %theta aggiornati
    end
    
    %y=y*(1/2); %scalatura
    
    %crea video
    
    h1=subplot(3,2,1);
    h12=plot(h1,t,x);
    h13=title('Somma dei segnali rumorosa (osservazioni)');
    h5=subplot(3,2,2);
    temp_so=rickerw(t,f,the(1,1),the(1,2))+rickerw(t,f,the(2,1),the(2,2))+rickerw(t,f,the(3,1),the(3,2));
    h52=plot(h5,t,temp_so);
    h53=title(strcat('Somma dei segnali stimati all''iterazione',num2str(p)));
    
    h2=subplot(3,2,3);
    h22=plot(h2,t,y(:,1));
    h23=title(strcat('Stima segnale 1 all''iterazione',num2str(p)));
    h3=subplot(3,2,4);
    h32=plot(h3,t,y(:,2));
    h33=title(strcat('Stima segnale 2 all''iterazione',num2str(p)));
    h4=subplot(3,2,5);
    h42=plot(h4,t,y(:,3));
    h43=title(strcat('Stima segnale 3 all''iterazione',num2str(p)));
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
    
    
    %M-step
    %ho idea che serva una ricerca esaustiva
    
    for k=1:3
        %costruzione della funzione di verosimiglianza, da minimizzare
        for i=1:length(A)
            for l=1:length(tau)
                cost(i,l,k)=cost(i,l,k)+(1/sign)*((y(:,k)-rickerw(t,f,A(i),tau(l))')'*(y(:,k)-rickerw(t,f,A(i),tau(l))'));
            end
        end
        %cerco il minimo
        [min_val ind_riga]=min(cost(:,:,k)); %ind_riga è un vettore, contiene la riga col minimo
        %per ogni colonna
        [min_v ind_col]=min(min_val);
        Amin=A(ind_riga(ind_col));
        tau_m=tau(ind_col);
        the(k,:)=[Amin tau_m];
        %the(:,1)=the(:,1)*(1/2); %scalatura
    end
    
%     new=imagesc(tau,A,cost(:,:,3));
%     new2=colorbar;
%     new21=xlabel('Ritardo [s]');
%     new22=ylabel('Ampiezza');
%     new23=title(strcat('Funzione costo all''iterazione ',num2str(p),'per il segnale 1'));
%     F2 = getframe(fig2);
%     aviobj2 = addframe(aviobj2,F2);
    
end

close(fig);
close(fig2);

aviobj = close(aviobj);
aviobj2 = close(aviobj2);
    
%end

%La selezione nella costruzione del segnale somma fa si che a convergenza
%il segnale separato venga sommato a sè stesso, risultando dunque con
%ampiezza raddoppiata; per riportare le cose alla normalità lo si scala

% y=y*(1/2);
% the(:,1)=the(:,1)*(1/2);

%plotting
% figure,subplot(2,2,1),plot(t,x),title('Segnale ricevuto rumoroso'),...
%     subplot(2,2,2),plot(t,y(:,1)),title('Stima segnale 1'),...
%     subplot(2,2,3),plot(t,y(:,2)),title('Stima segnale 2'),...
%     subplot(2,2,4),plot(t,y(:,3)),title('Stima segnale 3'),...

%Video dell'analisi via STFT, come varia lo spettro nel tempo

 
% 
% for k=1:size(TRA,2)
%     po=plot(fk,abs(TRA(:,k)));
%     h1=xlabel('Frequenza [Hz]');
%     h2=ylabel('Ampiezza');
%     h3=title('Evoluzione dello spettro della Ricker nel tempo');
%     F = getframe(fig);
%     aviobj = addframe(aviobj,F);
% end
% 
% close(fig);
% aviobj = close(aviobj);
    

