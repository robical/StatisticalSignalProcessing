%funzione che calcola la DFT in forma matriciale

% x= sequenza in ingresso
% m= numero di campioni dft da calcolare

function [X,fk]=dft(x,m,fc)

n=length(x);

%zero padding
if(n<m)
    x=[x zeros(1,m-n)];
end

W=zeros(m); %matrice DFT
for k=1:m
    for i=1:m
        W(k,i)=exp(-j*((2*pi)/m)*k*i);
    end
end

X=x*W.';

if(mod(m,2)==0)
    fk=-(m/2-1)*fc/m:fc/m:(m/2)*fc/m;
    X=[X(m/2+2:m) X(1:m/2+1)];
else
    fk=-((m-1)/2)*fc/m:fc/m:((m-1)/2)*fc/m;
    X=[X(((m-1)/2)+2:m) X(1:((m-1)/2)+1)];
end



        