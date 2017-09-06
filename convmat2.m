%Modello filtraggio lineare: identificazione

%crea la matrice di convoluzione con il segnale in ingresso
% u = segnale in ingresso (intero)
% h = filtro
% n = cardinalità dell'insieme delle osservazioni
% n0 = che campione sto considerando come campione zero?


function U=convmat2(u,h,n0,n)


p=length(h);




U=zeros(n,p);

ind_min=0;
ind_max=n-1;
for k=1:p
    U(:,k)=u(n0+ind_min:1:n0+ind_max)';
    ind_min=ind_min-1;
    ind_max=ind_max-1;
end
