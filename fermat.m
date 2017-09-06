%Funzione che calcola i percorsi intermedi BS-superficie e superficie-MS
%con raggio di Fermat --> l1 ed l2 rispettivamente
% d = distanza in piano BS-MS
% hBS =altezza BS in metri
% hMS = altezza MS in metri



function [l1,l2]=fermat(hBS,hMS,d)

    df1=((-d*(hBS^2))+(hMS*hBS*d))/(hMS^2-hBS^2);
    df2=((-d*(hBS^2))-(hMS*hBS*d))/(hMS^2-hBS^2);
    
    if (df1>0)
        df=df1;
    elseif (df2>0)
        df=df2;
    else
        df=min(df1,df2);
    end
    
    l1=sqrt(hBS^2 + df^2);
    l2=sqrt((d-df)^2 + hMS^2);