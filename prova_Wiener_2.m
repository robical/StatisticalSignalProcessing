%% Deconvoluzione MMSE con filtro di Wiener

%Devo conoscere:

% 1) filtro lineare attraversato
% 2) autocorrelazione delle osservazioni
% 3) traformata Z della sequenza di ingresso, almeno approssimativamente

clear all;
close all;
clc;


% Rumore bianco in ingresso

N=100;

%sig_w=2; %potenza del rumore ingresso

%w=sig_w*randn(N,1); %vettore colonna


n=0:N-1;

%f=10;

sig_w=10;

w=sig_w*randn(N,1); %segnale, vettore colonna

%Filtro lineare attraverso il quale passa il rumore
h=[1 -1 1 -1 1]'; %5 campioni

M=length(h);



%rumore bianco additivo

G=N+M-1;

sig_z=0.1; %varianza del rumore additivo

z=sig_z*randn(G,1);



%osservazioni

x=conv(w,h)+z;

X=fftshift(fft(x,G));

H=fftshift(fft(h,G));

%Conosco la DSP del processo di segnale, che è a banda piatta e vale sig_w

%stimo numeratore del filtro di Wiener

num=(sig_w.*conj(H));

den=(sig_w.*(H.*conj(H)))+sig_z;

A_stim=num./den;

asse=(-G/2+1)*((2*pi)/G):((2*pi)/G):((2*pi)/G)*G/2;

figure,subplot(1,2,1),plot(asse,abs(A_stim),'k'),title('Modulo del filtro di Wiener'),...
    subplot(1,2,2),plot(asse,angle(A_stim),'r'),title('Fase del filtro di Wiener')

figure,subplot(1,2,1),plot(asse,abs(W),'k'),title('Modulo del rumore in ingresso'),...
    subplot(1,2,2),plot(asse,angle(W),'r'),title('Fase del rumore in ingresso')

W_stim=X.*A_stim;

figure,subplot(1,2,1),plot(asse,abs(W_stim),'k'),title('Modulo del segnale in ingresso stimato'),...
    subplot(1,2,2),plot(asse,angle(W_stim),'r'),title('Fase del segnale in ingresso stimato')