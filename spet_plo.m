% funzione che plotta lo spettro in modulo e fase rappresentato dai valori
% della trasformata Z di un filtro con zeri e poli rappresentati dai
% polinomi da passare in ingresso Az(zeri) Bp(poli), calcolata sul cerchio
% unitario




function [Spettro,f]=spet_plo(Az,Bp,fc)

%Potrei a questo punto valutare la fdt per frequenze continue, z=exp(j*2*pi*(f/fc))

%Traccio la funzione di trasferimento corrispondente ad un filtro con
%trasformata Z che ha questi zeri



M=2^11; %numero di pulsazioni normalizzate nelle quali valuto la FDT
f=-fc/2+fc/M:fc/M:fc/2;

for k=1:length(f)
A=1;
for i=1:length(Az)
    A=A*(1-Az(i)*exp(-1i*2*pi*(f(k)./fc)));
end

Num(k)=A;

B=1;
for i=1:length(Bp)
    B=B*(1-Bp(i)*exp(-1i*2*pi*(f(k)./fc)));
end

Den(k)=B;

clear A;
clear B;
end

Spettro=Num./Den;