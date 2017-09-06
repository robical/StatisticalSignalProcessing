%% Stimatore MMSE per tx binaria antipodale

clear all;
close all;
clc;



% 
% for i=1:1000 %numero osservazioni
% sigmaw=1;
% w=sigmaw*randn;
% the=randsrc;
% 
% x=the+w;
% 
% c(i)=x;
% end
% 
% %figure,hist(c,500),title('Istogramma delle osservazioni')
% 
% 
% log_px=log(1/sqrt(8*pi*sigmaw))-(c.^2+1)/sigmaw;
% 
% the_mmse=c(1:end-1)+sigmaw*(diff(log_px));

%figure,hist(the_mmse,500),title('Istogramma della stima di Theta MMSE')

%prova distribuzione delle osservazioni

x=-7:.01:7;

sigmaw=[0.1 0.2 0.5 1 3];
col=['r' 'k' 'b' 'm' 'g'];

figure(1),
for i=1:5

px=1/sqrt(8*pi*sigmaw(i))*(exp(-((x+1).^2/(2*sigmaw(i))))+exp(-((x-1).^2/(2*sigmaw(i)))));

logpx=log(px);

%stima della derivata

px_dx=diff(logpx)./0.01;

%stima del the_mmse

the_mmse=x(1:end-1)+sigmaw(i)*px_dx;

w_mmse(i,:)=x(1:end-1)-the_mmse;
% figure,subplot(1,2,1),plot(x,px,'r'),title('Pdf delle osservazioni (f chiusa)'),...
%     subplot(1,2,2),plot(x,logpx,'k'),title('Log pdf osservazioni')

%figure,subplot(1,2,1),plot(x(1:end-1),px_dx,'b'),title('Stima della derivata logaritmica di Px'),...
    
    hold on,plot(x(1:end-1),the_mmse,col(i)),axis([min(x)-1 max(x)+1 min(x)-1 max(x)+1]),...
    title('Stimatore MMSE'),xlabel('x'),ylabel('E[the|x]'),hold on,
end

hold on, plot(x(1:end-1),x(1:end-1),'k-.'),legend('s_w=0.1','s_w=0.3','s_w=0.5','s_w=1','s_w=3','lineare')

%stimatore MMSE del rumore
figure,
for i=1:5
    hold on,plot(x(1:end-1),w_mmse(i,:),col(i)),axis([min(x)-1 max(x)+1 min(x)-1 max(x)+1]),...
    title('Stimatore MMSE rumore'),legend('s_w=0.1','s_w=0.3','s_w=0.5','s_w=1','s_w=3'),xlabel('x'),ylabel('E[w|x]'),hold on,
end
hold on, plot(x(1:end-1),x(1:end-1),'k-.')

%% Tx multilivello e stima MMSE
clear all;
close all;
clc;


%in questo caso il segnale è multilivello, ad esempio 4 livelli: -3 -1 1 3

x=-7:.01:7;

sigmaw=[0.1 0.2 0.5 1 3];
col=['r' 'k' 'b' 'm' 'g'];

%figure(1),
%for i=1:5
i=3;
px=1/sqrt(32*pi*sigmaw(i))*(exp(-((x+1).^2/(2*sigmaw(i))))+exp(-((x-1).^2/(2*sigmaw(i))))+exp(-((x+3).^2/(2*sigmaw(i))))+exp(-((x-3).^2/(2*sigmaw(i)))));

logpx=log(px);

%stima della derivata

px_dx=diff(logpx)./0.01;

%stima del the_mmse

the_mmse=x(1:end-1)+sigmaw(i)*px_dx;

w_mmse(i,:)=x(1:end-1)-the_mmse;
figure,subplot(1,2,1),plot(x,px,'r'),title('Pdf delle osservazioni (f chiusa)'),...
     subplot(1,2,2),plot(x,logpx,'k'),title('Log pdf osservazioni')

figure,subplot(1,2,1),plot(x(1:end-1),px_dx,'b'),title('Stima della derivata logaritmica di Px'),...
    
    figure,hold on,plot(x(1:end-1),the_mmse,col(i)),axis([min(x)-1 max(x)+1 min(x)-1 max(x)+1]),...
    title('Stimatore MMSE'),xlabel('x'),ylabel('E[the|x]'),hold on,
%end

hold on, plot(x(1:end-1),x(1:end-1),'k-.'),legend('s_w=0.1','s_w=0.3','s_w=0.5','s_w=1','s_w=3','lineare')

%stimatore MMSE del rumore
figure,
%for i=1:5
    hold on,plot(x(1:end-1),w_mmse(i,:),col(i)),axis([min(x)-1 max(x)+1 min(x)-1 max(x)+1]),...
    title('Stimatore MMSE rumore'),legend('s_w=0.1','s_w=0.3','s_w=0.5','s_w=1','s_w=3'),xlabel('x'),ylabel('E[w|x]'),hold on,
%end
hold on, plot(x(1:end-1),x(1:end-1),'k-.')

%% Mistura Gaussiana (GMM)
clear all;
close all;
clc;

%Definisco le 2 pdf che andranno a comporre la mistura, e la prob
%dell'evento di scelta

p=0.5;

mu_a=3; %media
sig_a=5; %varianza

mu_b=1;
sig_b=2;

sig_w=0.5;

%Simulazione di mistura gaussiana

binz=-20:0.1:20;

N=zeros(size(binz));

Nu=10000;

for k=1:200

for i=1:Nu

res=binornd(1,p);

if (res==0)
    a=mu_a+sig_a*randn; %estrazione da vc alfa
    x(i)=a+sig_w*randn;
else
    b=mu_b+sig_b*randn; %estrazione da vc beta
    x(i)=b+sig_w*randn;
end

end

N=N+histc(x,binz);

end

N=N./sum(N); %stima della pdf delle osservazioni

figure,plot(binz,N),title('Stima della pdf delle osservazioni vs pdf calcolata'),hold on,


%Calcolo teorico delle pdf alfa e beta

clear x;

x=binz;

Ga_x=(1/(sqrt(2*pi*sig_a)))*exp(-((x-mu_a).^2)/(2*sig_a));

Gb_x=(1/(sqrt(2*pi*sig_b)))*exp(-((x-mu_b).^2)/(2*sig_b));

%Definisco la pdf del rumore (gaussiano bianco, dunque a valore medio nullo e incorrelato)


Gw_x=(1/(sqrt(2*pi*sig_w)))*exp(-(x.^2)/(2*sig_w));

%Calcolo la pdf teorica delle osservazioni con il GMM

p=0.5;

pxa_rid=conv(Ga_x,Gw_x,'same'); %va rinormalizzata

pxa_rid=pxa_rid./sum(pxa_rid);

pxb_rid=conv(Gb_x,Gw_x,'same'); %va rinormalizzata

pxb_rid=pxb_rid./sum(pxb_rid);

%previsione teorica della pdf delle osservazioni, da confrontare poi con la
%simulazione

px=p*pxa_rid+(1-p)*pxb_rid; %va rinormalizzata

px=px./sum(px);

%La plotto sopra, e vedo se coincide con la simulazione

plot(binz,px,'r'),legend('Blu - pdf stimata','Rosso - pdf teorica'),ylabel('P_x(x)'),...
    xlabel('x')


