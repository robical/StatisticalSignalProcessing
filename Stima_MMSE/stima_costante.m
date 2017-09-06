%% Stima MMSE per processi gaussiani e modello lineare (stima di costante)

%Nella stima bayesiana è fondamentale la scelta ragionata dell'utilizzo di
%una certa pdf a PRIORI, perchè è su di essa che si basa poi la stima del
%parametro in caso di rumore molto alto e/o poche osservazioni

clear all;
close all;
clc;


%Modello lineare

N=10; %numero di osservazioni
P=3; %valori assunti dal parametro durante TUTTA l'osservazione del fenomeno

sig_w=3; %varianza (potenza) di rumore

%Definisco il parametro come una vc gaussiana con un certo valore medio e
%varianza - acquisisco almeno 10 osservazioni per ogni valore del parametro

med_par=3;
sig_par=1;

for k=1:100

for i=1:P


w=sig_w*randn(N,1); %a media nulla

the=sig_par*(randn(P,1)+med_par);

x(:,i)=the(i)+w; %ho 10 osservazioni per ogni valore del parametro, cioè per ogni 
%colonna della matrice osservazioni

end

%Stimatore MMSE del valore della costante
alfa=(sig_w/N)/(sig_w/N+sig_par);

for i=1:P
the_mmse(k,i) = alfa*med_par+(1-alfa)*((ones(1,N)*x(:,i))/N);
end
end

the_mmse_med=mean(the_mmse');

figure,plot(1:100,the_mmse_med) %mostro la media delle 3 stime MMSE effettuate
%su 3 diverse realizzazioni del parametro, ricordando che ho solo 10
%osservazioni per ogni realizzazione del parametro

media_stima=mean(the_mmse_med);

disp(strcat('La stima del parametro è: ',num2str(media_stima)))