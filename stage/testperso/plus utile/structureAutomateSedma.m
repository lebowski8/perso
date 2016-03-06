% Test structure d'automate
clear all
%close all
home

Ts=0.005;                % Pas de temps
thetaK=18;                  % Temperature initiale
Tmax=22;              % Température à ne pas dépasser en maxi
Tmin = 15;             % Température à ne pas dépasser en mini
tsim=5;                % Temps de simulation

t = 0:Ts:tsim-Ts;       % Vecteur de temps

in=(Tmax+1)*ones(tsim/Ts,1);  % Entree positive
X=zeros(tsim/Ts,1);           % Vecteur d'état
theta=zeros(tsim/Ts,1);           % Vecteur de sortie
X(1)=thetaK;                      % Etat initial (pris comme la sortie initiale)


%% Gestion automate

%Declaration des etats et initialisation
%unsigned int ES,EP;	% Variables Etat suivant et Etat present
ES=0;
EP=0;

% Correspondance entre les nom d'etats et le vecteur entier : 
%  E0 : 0
%  E1 : 1
%

%% Debut du bloc F de la realisation

for k=1:tsim/Ts
    switch(EP) 
        % code pour Ep = E0
        case 0 
            if(thetaK >= Tmax)
                ES = 1;
            end
        % code pour Ep = E1
        case 1 
            if(thetaK <= Tmin)		
                ES = 0;		
            end
    end

    EP=ES;
    %Fin du bloc F de la realisation

    %Debut de realisation du bloc G : 

    if EP == 0 
       A=-1;B=1;C=1;D=0;       % Représentation d'état (continue)
%         A = [0 1 ; -1 -2];
%         B = [0 ; 0];
%         C = [1 0];
%         D = [0];
        
        
        
       sys1=ss(A,B,C,D);
       sys1=c2d(sys1,Ts);      % Passage en discret
       X(k+1)=sys1.A*X(k)+sys1.B*in(k); 
       theta(k)=sys1.C*X(k)+sys1.D*in(k);
       thetaK=theta(k);
    else
        marche = 0;
    end

    if EP == 1 
        A=-1;B=0;C=1;D=0;        % Système sans entrée (switch si T>Tlim)
        sys2=ss(A,B,C,D);
        sys2=c2d(sys2,Ts);
        X(k+1)=sys2.A*X(k)+sys2.B*in(k);
        theta(k)=sys2.C*X(k)+sys2.D*in(k);
        thetaK=theta(k);
    else
        arret = 0;
    end
    %Fin de realisation du bloc G 

end




%% Affichage 
figure
plot(t',theta)
grid on
xlabel('temps')
ylabel('temperature')


