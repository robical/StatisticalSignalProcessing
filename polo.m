function [y]=polo(x,z0)


%OK se polo � a fase minima abs(z0)<1
if(abs(z0)<1)
    
    in=1e-5; %valore che in realt� � di y(-1), serve come impulso ad alimentare
    %il meccanismo ricorsivo che alimenta il polo
    y(1)=x(1)+z0*in; %impuso in ingresso che alimenta il polo (cio� X(z) a banda piatta, 1)
    y(2)=x(2)+z0*y(1);
    
    i=2;
    while(i<length(x))
        y(i+1)=x(i+1)+z0*y(i);
        i=i+1;
    end
    
else
    in=1e-5; %stavolta � il campione che alimenta il circuito ricorsivo n=1
    %trover� la sequenza convergente anticausale, cio� che va da 0 a -n
    y(1)=(1/z0)*(in-x(1));
    y(2)=(1/z0)*(y(1)-x(2));

    i=2;
    while(i<length(x))
        y(i+1)=(1/z0)*(y(i)-x(i));
        i=i+1;
    end
    y=fliplr(y);
end



