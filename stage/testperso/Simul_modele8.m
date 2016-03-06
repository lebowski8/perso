%% Cplex + MatLab
clear all
close all
home

load('modele.mat')

nEtat = 11; % nombre d'état du système
nEntre = 2; % nombre d'entrée
nSortie = 3; % nombre de sortie

Ts=0.5;                % Pas de temps en seconde
% Ts=2.5;                % Pas de temps en seconde
tsim=6*Ts;                % Temps de simulation
t = 0:Ts:tsim-Ts;       % Vecteur de temps
i=1;

eps = 50; % variation d'altitude autour du mode maintien
h(1) = 10000*0; % altitude initiale
Va = 235*0; % vitesse initiale
X(:,1) = zeros(nEtat,1);
X(4,1) = Va;
X(8,1) = h(1);




% u = [ 500*ones(1,length(t)); 300*ones(1,length(t)/2)
% 300*ones(1,length(t)/2)]; % ce plan est instable pour le mode maintien :
% on est pas dans l'ellipse formée par le plan (Va,alpha)
% u = [ 235*ones(1,length(t)); linspace(200,300,length(t)/2) 300*ones(1,length(t)/2)];
u = [ 20*ones(1,length(t)); 100*ones(1,length(t))];

% u = [235*ones(1,length(t)) ; 1499.9 1499.9 1499.9 1126.9 1499.9]; % test planif

% u = [0 235.1 350 350 350 0 ; 50 50.787 70.41 23.775 105.51 200];

tswitch = -1;

%% On fait en sorte que tous les modèles soit de même taille : nEtat = 11 or certains modèle n'ont que 10 état
        sys_descent11 = ss();
        sys_descent11.A = [sys_descent.A zeros(10,1) ; zeros(1,11)];
        sys_descent11.B = [sys_descent.B ; zeros(1,2)];
        sys_descent11.C = [sys_descent.C zeros(3,1)];
        sys_descent11.D = sys_descent.D;
        sys_descent11.Ts = sys_descent.Ts;
       % clear sys_descent;
        sys_descent = sys_descent11;
        clear sys_descent11;
        
        sys_climb11 = ss();
        sys_climb11.A = [sys_climb.A zeros(10,1) ; zeros(1,11)];
        sys_climb11.B = [sys_climb.B ; zeros(1,2)];
        sys_climb11.C = [sys_climb.C zeros(3,1)];
        sys_climb11.D = sys_climb.D;
        sys_climb11.Ts = sys_climb.Ts;
        %clear sys_climb;
        sys_climb = sys_climb11;
        clear sys_climb11;
        

%% création du système échantillonnée : Xk+n = f(Xk)
n = 5; % passage de 0.5 s à n*0.5s d'échantillonage

% sys_hold
sys_holdk5 = sys_hold;
sys_holdk5.Ts = sys_hold.Ts*n;

% Xk+1 = ...
sys_holdk5.a = sys_hold.a^n;
b = sys_hold.a^0;
for k = 1:n-1
    b = b + sys_hold.a^k;
end
sys_holdk5.b = b*sys_hold.b;

% Yk+1 = ...
sys_holdk5.c = sys_hold.c*sys_hold.a^n;
d = sys_hold.a^0;
for k = 1:n-1
    d = d + sys_hold.a^k;
end
d = d*sys_hold.b;
d = sys_hold.c*d;
d = d + sys_hold.d;

    % sys_climb
    sys_climbk5 = sys_climb;
    sys_climbk5.Ts = sys_climb.Ts*n;


    % Xk+1 = ...
    sys_climbk5.a = sys_climb.a^n;
    b = sys_climb.a^0;
    for k = 1:n-1
        b = b + sys_climb.a^k;
    end
    sys_climbk5.b = b*sys_climb.b;

    % Yk+1 = ...
    sys_climbk5.c = sys_climb.c*sys_climb.a^n;
    d = sys_climb.a^0;
    for k = 1:n-1
        d = d + sys_climb.a^k;
    end
    d = d*sys_climb.b;
    d = sys_climb.c*d;
    d = d + sys_climb.d;
    
% sys_descent
sys_descentk5 = sys_descent;
sys_descentk5.Ts = sys_descent.Ts*n;


% Xk+1 = ...
sys_descentk5.a = sys_descent.a^n;
b = sys_descent.a^0;
for k = 1:n-1
    b = b + sys_descent.a^k;
end
sys_descentk5.b = b*sys_descent.b;

% Yk+1 = ...
sys_descentk5.c = sys_descent.c*sys_descent.a^n;
d = sys_descent.a^0;
for k = 1:n-1
    d = d + sys_descent.a^k;
end
d = d*sys_descent.b;
d = sys_descent.c*d;
d = d + sys_descent.d;

%% Définition des systèmes et étude de stabilité : 
% figure
% hold on
% grid on

% HOLD EP = 0 :
    [P0 x0 y0 P_hold] = domaineLyap(sys_hold,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD
% CLIMB EP = 1 : 
    [P1 x1 y1 P_climb] = domaineLyap(sys_climb,'discret');
% DESCENT EP = -1 : 
    [P2 x2 y2 P_descent] = domaineLyap(sys_descent,'discret');
    
%% Etude de la stabilité des système k5

% HOLD EP = 0 :
    [P0 x0 y0 P_holdk5] = domaineLyap(sys_holdk5,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD
%     [P0k10 x0k10 y0k10] = domaineLyap(sys_holdk10,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD

% CLIMB EP = 1 : 
    [P1 x1 y1 P_climbk5] = domaineLyap(sys_climbk5,'discret');
% DESCENT EP = -1 : 
    [P2 x2 y2 P_descentk5] = domaineLyap(sys_descentk5,'discret');
    
%% Génération des Systèmes pour le format de donnée de Cplex
clc    
    
    mode = {-1, 0, 1};
    sysName = {'_Descent', '_Hold', '_Climb'};
    modes = {sys_descent, sys_hold, sys_climb};
    
    
    disp(['m = ', int2str(nEntre),';'])
    disp(['r = ', int2str(nSortie),';'])
    disp(['n = ', int2str(nEtat),';'])
    
    
    disp('// Définition des systèmes')
    
    for k = 1:length(modes)
        A = modes{k}.A;
        B = modes{k}.B;
        C = modes{k}.C;
        D = modes{k}.D;
        nEtat = length(A);
        mat = {A,B,C,D};
        matName = {'A','B','C','D'};
        for m = 1:length(mat)
            disp([matName{m}, sysName{k}, ' = ['])
            M = mat{m};
            for i = 1:size(M,1)
                chaine = '[';
                for j = 1:(size(M,2)-1)
                    chaine = [chaine , num2str(M(i,j),10),' , '];
                end
                disp(['         ',chaine, num2str(M(i,size(M,2)),10), '],'])
            end
            disp('];')
        end
    end
    

%% Gestion de l'automate

%Declaration des etats et initialisation
u;

EP = 0; % Etat present
ES = EP; % Etat suivant

% sys_hold = sys_holdk5;
% sys_climb = sys_climbk5;
% sys_descent = sys_descentk5;


for k=1:length(t)
    Vac = u(1,k);
    hc = u(2,k);
    % Debut du bloc F de la realisation
    switch(EP) 
        % code pour Ep = -1 : descent
        case -1             
            if((hc-h(k)) > -eps)
                disp('Switch de DESCENTE vers MAINTIEN')
                temps = t(k)
                ES = 0;
                tswitch(i) = t(k);
                affSwitch{i} = 'DESCENTE -> MAINTIEN';
                i= i+1;
            elseif ((hc-h(k)) >= eps)
                disp('Switch de DESCENTE vers MONTEE')
                temps = t(k)
                ES = 1;
                tswitch(i) = t(k);
                affSwitch{i} = 'DESCENTE -> MONTEE';
                i= i+1;
            else
                disp('Reste en DESCENTE')
                temps = t(k)
                ES = -1;
                %tswitch(i) = t(k);
                %i= i+1;
            end
        % code pour Ep = 0 : maintien
        case 0             
            if((hc-h(k))<=-eps)
                disp('Switch de MAINTIEN vers DESCENTE')
                temps = t(k)
                ES = -1;
                tswitch(i) = t(k);
                affSwitch{i} = 'MAINTIEN -> DESCENTE';
                i= i+1;
            elseif ((hc-h(k))>=eps)
                disp('Switch de MAINTIEN vers MONTEE')
                temps = t(k)
                ES = 1;
                tswitch(i) = t(k);
                affSwitch{i} = 'MAINTIEN -> MONTEE';
                i= i+1;
            else
                disp('Reste en MAINTIEN')
                temps = t(k)
                ES = 0;
                %tswitch(i) = t(k);
                %i= i+1;
            end
        % code pour Ep = 1 : monté
        case 1 
            if((hc-h(k)) < eps)
                disp('Switch de MONTEE vers MAINTIEN')
                temps = t(k)
                ES = 0;	
                tswitch(i) = t(k);
                affSwitch{i} = 'MONTEE -> MAINTIEN';
                i= i+1;
            elseif ((hc-h(k)) <= -eps)
                disp('Switch de MONTEE vers DESCENTE')
                temps = t(k)
                ES = -1;	
                tswitch(i) = t(k);
                affSwitch{i} = 'MONTEE -> DESCENTE';
                i= i+1;
            else
                disp('Reste en MONTEE')
                temps = t(k)
                ES = 1;
                %tswitch(i) = t(k);
                %i= i+1;
            end
    end

    EP=ES
    %Fin du bloc F de la realisation

    %Debut de realisation du bloc G : 
    if EP == -1 % sys_descent     
        A = [sys_descent.A zeros(10,1) ; zeros(1,11)];
        B = [sys_descent.B ; zeros(1,2)];
        C = [sys_descent.C zeros(3,1)];
        D = sys_descent.D;
       X(:,k+1)=A*X(:,k)+B*u(:,k); 
       Y(:,k)=C*X(:,k)+D*u(:,k);
        
%         disp('MODE DESCENT')
%         disp('stable à Xk ?')
%         inEllipse(P_descent,X(:,k))
%         disp('stable à Xk+1 ?')
%         inEllipse(P_descent,X(:,k+1))
        
%        X(:,k+1)=sys_descent.A*X(:,k)+sys_descent.B*u(:,k); 
%        Y(:,k)=sys_descent.C*X(:,k)+sys_descent.D*u(:,k);
       Yk=Y(:,k);
        Va(k) = Yk(1);
        h(k+1) = Yk(2);
        alpha(k) = Yk(3);
    end
    
    if EP == 0 % sys_hold
%         A = sys_hold.A;
%         B = sys_hold.B;
%         C = sys_hold.C;
%         D = sys_hold.D;
       X(:,k+1)=sys_hold.A*X(:,k)+sys_hold.B*u(:,k); 
       Y(:,k)=sys_hold.C*X(:,k)+sys_hold.D*u(:,k);
%        disp('MODE HOLD')
%         disp('stable à Xk ?')
%         inEllipse(P_hold,X(:,k))
%         disp('stable à Xk+1 ?')
%         inEllipse(P_hold,X(:,k+1))
       Yk=Y(:,k);
        Va(k) = Yk(1);
        h(k+1) = Yk(2);
        alpha(k) = Yk(3);
    end

    if EP == 1 % sys_climb
        A = [sys_climb.A zeros(10,1) ; zeros(1,11)];
        B = [sys_climb.B ; zeros(1,2)];
        C = [sys_climb.C zeros(3,1)];
        D = sys_climb.D;
       X(:,k+1)=A*X(:,k)+B*u(:,k); 
       Y(:,k)=C*X(:,k)+D*u(:,k);
       
%        disp('MODE CLIMB')
%         disp('stable à Xk ?')
%         inEllipse(P_climb,X(:,k))
%         disp('stable à Xk+1 ?')
%         inEllipse(P_climb,X(:,k+1))       
%         X(:,k+1)=sys_climb.A*X(:,k)+sys_climb.B*u(:,k); 
%         Y(:,k)=sys_climb.C*X(:,k)+sys_climb.D*u(:,k);
        Yk=Y(:,k);
        Va(k) = Yk(1);
        h(k+1) = Yk(2);
        alpha(k) = Yk(3);
    end
    %Fin de realisation du bloc G 
end


%% Affichage 

figure
plot(t,Va)
title('Va(t)')
grid on
hold on

if tswitch(1) == -1 
    disp('Aucun switch')
else
    for u = 1:length(tswitch)
        disp(['t = ' ,num2str(tswitch(u)), ' : ', affSwitch{u}])
        x = tswitch(u);
        y = min(Va):1:max(Va);
        plot(x,y,'r'); % affiche les instants de commutations
    end
end

figure
plot(t,h(1:length(t)))
title('h(t)')
grid on
hold on

if tswitch(1) == -1 
    disp('Aucun switch')
else
    for u = 1:length(tswitch)
        x = tswitch(u);
        y = min(h):1:max(h);
        plot(x,y,'r'); % affiche les instants de commutations
    end
end

figure
plot(t,alpha*180/pi)
title('alpha(t)')
grid on
hold on

if tswitch(1) == -1 
    disp('Aucun switch')
else
    for u = 1:length(tswitch)
        x = tswitch(u);
        y = min(alpha*180/pi):0.1:max(alpha*180/pi);
        plot(x,y,'r'); % affiche les instants de commutations
    end
end

%%
hcInit = 1e4;
hc = hcInit + 500;
sim('Discrete_Closed_Loop_Model_74_prelude_simu');
Va(end)