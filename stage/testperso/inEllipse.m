function ok = inEllipse(P,x)
%% inEllipse permet de savoir si le point x est dans l'ellipse définit par P 
x = x(1:length(P));
    if x'*P*x <= 1
        ok = true;
    else
        ok = false;
    end

end