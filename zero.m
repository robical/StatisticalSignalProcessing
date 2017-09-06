%funzione che genera la sequenza complessa pertinente ad una certa
%trasformata Z contenente solo zeri


function [A]=zero(z)

A=[1 -z(1)];
for i=2:length(z)
    A=conv(A,[1 -z(i)]);
end