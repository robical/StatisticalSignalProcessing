function [H]=Ztra(h,N)

syms z;

M=length(h);

h=[h' zeros(1,N-M)]';

for i=0:N-1
    H=h(i+1)*z^-i; %trasf Z del filtro
end
