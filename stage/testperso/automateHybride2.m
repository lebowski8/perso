% Test structure d'automate
clear all
close all
home


Ts=0.05;                % Pas de temps
Tmax=22;              % Température à ne pas dépasser en maxi
Tmin = 15;             % Température à ne pas dépasser en mini
tsim=5;                % Temps de simulation
t = 0:Ts:tsim-Ts;       % Vecteur de temps
i=1;

u=(Tmax+1)*ones(tsim/Ts,1);   % Entree positive
X=zeros(2,tsim/Ts);           % Vecteur d'état
Y=zeros(2,tsim/Ts);           % Vecteur de sortie

Yk = [0.5 ; 0.5];
% Yk = [1e3 ; 5e2];
X(:,1) = Yk;                  % Etat initial (pris comme la sortie initiale)
tswitch = -1;


%% Définition des systèmes : 
% E0 :
    A = [0 1 ; -1 -2];
    B = [0 ; 0];
    C = [1 0 ; 0 1];
    D = [0 ; 0];

    sys0c=ss(A,B,C,D);
    [P0 x0 y0] = domaineLyap(sys0c);
    sys1=c2d(sys0c,Ts);      % Passage en discret
    %domaineLyap(sys1,'discret');
% E1 : 
    A = [0 5 ; -1 -2];
    B = [0 ; 0];
    C = [1 0 ; 0 1];
    D = [0 ; 0];

    sys1c=ss(A,B,C,D);
    sys2=c2d(sys1c,Ts);

    [P1 x1 y1] = domaineLyap(sys1c);
    %domaineLyap(sys2,'discret');
    
%% Domaines de stabilité global

figure
hold on
grid on
axis equal
title('Domaines de stabilités du système global')
xlabel('x1')
ylabel('x2')

plot(x0,y0)
plot(x1,y1,'r')

legend('E0', 'E1')

%% Gestion de l'automate

%Declaration des etats et initialisation

EP = 1; % Etat present
ES = EP; % Etat suivant

for k=1:tsim/Ts
    % Debut du bloc F de la realisation
    switch(EP) 
        % code pour Ep = E0
        case 0             
            if((Yk(1) <= 0.2) && (inEllipse(P1,Yk) == 1))
%                 if (inEllipse(P1,Yk) == 0) %permet de savoir si le système est dans l'ellipses du mode suivant
%                     disp('Yk non stable pour le mode suivant (E1)')
%                     disp(Yk)
%                     disp(t(k))
%                     %break;
%                 end
                disp('Switch de E0 vers E1')
                temps = t(k)
                ES = 1;
                tswitch(i) = t(k);
                i= i+1;
            end
        % code pour Ep = E1
        case 1 
            if((Yk(1) >= 0.6) && (inEllipse(P0,Yk) == 1))
%                 if (inEllipse(P0,Yk) == 0)
%                     disp('Yk non stable pour le mode suivant (E0)')
%                     disp(Yk)
%                     disp(t(k))
%                     %break;
%                 end
                disp('Switch de E1 vers E0')
                temps = t(k)
                ES = 0;	
                tswitch(i) = t(k);
                i= i+1;
            end
    end

    EP=ES;
    %Fin du bloc F de la realisation

    %Debut de realisation du bloc G : 

    if EP == 0 
        
       if (inEllipse(P0,Yk) == 0) % Permet de savoir si l'état est stable
          %disp('non stable');
          %break;
       end
       
       X(:,k+1)=sys1.A*X(:,k)+sys1.B*u(k); 
       Y(:,k)=sys1.C*X(:,k)+sys1.D*u(k);
       Yk=Y(:,k);

    end

    if EP == 1 
        
        if (inEllipse(P1,Yk) == 0)
          %disp('non stable'); 
          %break;
        end
        
        X(:,k+1)=sys2.A*X(:,k)+sys2.B*u(k);
        Y(:,k)=sys2.C*X(:,k)+sys2.D*u(k);
        Yk=Y(:,k);

    end
    %Fin de realisation du bloc G 
    plot(Yk(1),Yk(2),'g+')

end




%% Affichage 

figure
plot(t,Y)
grid on
xlabel('temps')
ylabel('état')
legend('x1','x2')
hold on

if tswitch(1) == -1 
    disp('Aucun switch')
else
    for u = 1:length(tswitch)
        x = tswitch(u);
        y = min(min(Y(2,:),Y(1,:))):0.01:max(max(Y(2,:),Y(1,:)));
%         y = min(min(Y(2,:)),min(Y(1,:))):1:max(max(Y(2,:)),max(Y(1,:)));
        plot(x,y,'r'); % affiche les instants de commutations
    end
end


