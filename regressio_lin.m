%regressione lineare

clear all;
close all;
clc;


fc=20;
%Intervallo non simmetrico
tmax=2;
tmin=0;

t=linspace(tmin,tmax,fc*(tmax-tmin));

a=5;
b=3;

y=a*t+b;

%figure,plot(t,y),title('Retta da stimare')

%Stima ML

%y=H*the+w

for k=1:30

sigma_w=k;

w=sigma_w*randn(size(y));

x(k,:)=y+w;

%figure,plot(t,x),title('Osservazioni con rumore')

H=zeros(length(t),2);

H=[ones(length(t),1) t'];

%covarianza del rumore è nota

Cw=sigma_w*eye(length(x(k,:))); %non conoscendo a priori le covarianze, le vc si suppongono incorrelate,
%commettendo un errore

the(:,k)=(H'*Cw^-1*H)^-1*H'*Cw^-1*x(k,:)';

%retta di regressione
%hold on, plot(t,the(1)+the(2)*t,'g')

%plotting della varianza della stima del parametro A
cova=(H'*Cw^-1*H)^-1;

var_a(k)=cova(1,1);
var_b(k)=cova(2,2);
end

figure,plot(10*log10(1:30),var_a,'r',10*log10(1:30),var_b,'b'),title('Varianza della stima dei parametri di una retta di regress, all''aumentare del rumore'),...
    xlabel('Potenza di rumore [dB]'),ylabel('Varianza del parametro'),legend('var. intercetta','var. coef. ango')
figure,plot(t,x(1,:),'*r'),hold on,plot(t,the(1,1)+the(2,1)*t,'g'),...
    title('Retta di regressione per potenza di rumore sigma=1'),xlabel('Tempo [s]'),ylabel('Ampiezza campioni')



fc=20;
%Intervallo simmetrico
tmax=1;
tmin=-1;

t=linspace(tmin,tmax,fc*(tmax-tmin));

a=5;
b=3;

y=a*t+b;

%figure,plot(t,y),title('Retta da stimare')

%Stima ML

%y=H*the+w

for k=1:30

sigma_w=k;

w=sigma_w*randn(size(y));

x=y+w;

%figure,plot(t,x),title('Osservazioni con rumore')

H=zeros(length(t),2);

H=[ones(length(t),1) t'];

%covarianza del rumore è nota

Cw=sigma_w*eye(length(x)); %non conoscendo a priori le covarianze, le vc si suppongono incorrelate,
%commettendo un errore

the=(H'*Cw^-1*H)^-1*H'*Cw^-1*x';

%retta di regressione
%hold on, plot(t,the(1)+the(2)*t,'g')

%plotting della varianza della stima del parametro A
cova=(H'*Cw^-1*H)^-1;

var_a_sim(k)=cova(1,1);
var_b_sim(k)=cova(2,2);
end

figure,plot(10*log10(1:30),10*log10(var_a),'r',10*log10(1:30),10*log10(var_a_sim),'b',10*log10(1:30),10*log10(var_a./var_a_sim),'k'),title('Confronto tra accuratezza nella stima dell''intercetta per misure simmetriche e non'),...
    xlabel('Potenza di rumore [dB]'),ylabel('Varianza della stima di a [dB]'),...
    legend('misure non simmetriche- accuratezza','misure simmetriche- accuratezza','guadagno')
figure,plot(10*log10(1:30),var_a_sim,'r',10*log10(1:30),var_b_sim,'b'),title('Varianza della stima dei parametri di una retta di regress, all''aumentare del rumore - SIM'),...
    xlabel('Potenza di rumore [dB]'),ylabel('Varianza del parametro'),legend('var. intercetta','var. coef. ango')