function [mod,fas,asse]=Zplot(X,M)

%Utilizzo:

% X= polinomio che rappresenta la trasformata Z della sequenza
% M=numero di campioni della DFT
% mod= modulo calcolato nei punti di campionamento equidistanti sul cerchio unitario
% fas= fase
% asse= asse delle pulsazioni normalizzate

g=1;
for i=-M/2+1:M/2 %pongo al centro lo zero
w=i*(2*pi)/M; %pulsazione ---> posizione angolare sul cerchio unitario sulla quale
%valuto
P=subs(X,z,exp(j*w));
modulo(g)=abs(P); %effettuo il campionamento regolare della trasformata 
fase(g)=angle(P);
g=g+1;
end

asse=(-M/2+1)*((2*pi)/M):((2*pi)/M):((2*pi)/M)*(M/2);
