% Test structure d'automate
clear all
close all
home

load('modele.mat')

nEtat = 11; % nombre d'état du système
nEntre = 2; % nombre d'entrée
nSortie = 3; % nombre de sortie

% Ts=0.5;                % Pas de temps en seconde
Ts=2.5;                % Pas de temps en seconde
tsim=7.5;                % Temps de simulation
t = 0:Ts:tsim-Ts;       % Vecteur de temps
i=1;

eps = 50; % variation d'altitude autour du mode maintien
h(1) = 1000; % altitude initiale
Va = 235; % vitesse initiale
X(:,1) = zeros(nEtat,1);
X(4,1) = Va;
X(8,1) = h(1);




% u = [ 500*ones(1,length(t)); 300*ones(1,length(t)/2)
% 300*ones(1,length(t)/2)]; % ce plan est instable pour le mode maintien :
% on est pas dans l'ellipse formée par le plan (Va,alpha)
% u = [ 235*ones(1,length(t)); linspace(200,300,length(t)/2) 300*ones(1,length(t)/2)];
% u = [ 235*ones(1,length(t)); 200*ones(1,length(t))];

u = [235*ones(1,length(t)) ; 0 0 1000]; % test planif

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
    [P0 x0 y0] = domaineLyap(sys_hold,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD
%     [P0k10 x0k10 y0k10] = domaineLyap(sys_holdk10,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD

% CLIMB EP = 1 : 
    [P1 x1 y1] = domaineLyap(sys_climb,'discret');
% DESCENT EP = -1 : 
    [P2 x2 y2] = domaineLyap(sys_descent,'discret');
    
%% Etude de la stabilité des système k10

% HOLD EP = 0 :
    [P0 x0 y0] = domaineLyap(sys_holdk5,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD
%     [P0k10 x0k10 y0k10] = domaineLyap(sys_holdk10,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD

% CLIMB EP = 1 : 
    [P1 x1 y1] = domaineLyap(sys_climbk5,'discret');
% DESCENT EP = -1 : 
    [P2 x2 y2] = domaineLyap(sys_descentk5,'discret');
    
%% Génération des Systèmes en CSP
clc
    m = 1;

    if m == -1
        sys = sys_descentk5;
        variableTemp = {'tempDescent'};
        nEtat = 10;
    elseif m == 0
        sys = sys_holdk5;
        nEtat = 11;
        variableTemp = {'tempHoldk5'};
    elseif m == 1
        sys = sys_climbk5;
        nEtat = 10;
        variableTemp = {'tempClimb'};    
    end
    
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    nTemp1 = ((nEtat+nEntre)*nEntre)*(nEtat);
    nTemp2 = ((nEtat+nEntre)*nEntre)*(nSortie);

    nTemp = nTemp1 + nTemp2;
    
    
    disp(['FloatVar[][] ',variableTemp{1},' = new FloatVar[step][',int2str(nTemp),'];'])
    disp(['for(int j = 0 ; j < step ; j++)'])
    disp(['    for(int i = 0 ; i < ', int2str(nTemp) ,' ; i++){'])
        disp(['        ',variableTemp{1},'[j][i] = new FloatVar(store, "',variableTemp{1},'"+i, -1e150, 1e150);']);
    disp('    }')
    
    disp(['for(int i = 1 ; i < step ; i++){'])

    for k = 1:nEtat % Calcul de Xk+1 = A*Xk + B*Uk
        disp('    // ##############################')

        disp('    // calcul des A(k,:)*Xk')
        for j = 1:nEtat % pour toutes les colonnes de A            
            %disp((k-1)*(nEtat+nEntre)*2 + j)
            ind = (k-1)*(nEtat+nEntre)*2 + j -1;
            disp(['    store.impose( new PmulCeqR(vars[i-1][',int2str(7+j),'],',num2str(A(k,j)),',',variableTemp{1},'[i-1][',int2str(ind),']) );'])
        end

        disp('    // calcul des B*Uk')
        for j = 1:nEntre           
            %disp(k*nEtat + (k-1)*(nEtat+nSortie+1) + j)
            ind = k*nEtat + (k-1)*(nEtat+nSortie+1) + j -1;
            disp(['    store.impose( new PmulCeqR(vars[i-1][',int2str(j-1),'],',num2str(B(k,j)),',',variableTemp{1},'[i-1][',int2str(ind),']) );'])            
        end

        disp('    // calcul des sommes (a11*x1k + a12*x2k ...)')
        ind = (k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + 1 -1;
        disp(['    store.impose( new PeqQ(',variableTemp{1},'[i-1][',int2str(ind),'],',variableTemp{1},'[i-1][',int2str((k-1)*(nEtat+nEntre)*2 + 1 -1),']) );'])            
        for j = 2:(nEtat+nEntre)            
            %disp((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + j) 
            ind = (k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + j -1;
            disp(['    store.impose( new PplusQeqR(',variableTemp{1},'[i-1][',int2str(ind-1),'],',variableTemp{1},'[i-1][',int2str((k-1)*(nEtat+nEntre)*2 + j -1),'],',variableTemp{1},'[i-1][',int2str(ind),']) );'])            
        end

%             disp(['    // vars[i][',int2str(7+k),'] = A(1,:)*x(k) + B(1,:)*u(k)'])
        disp(['    // m(i-1) = ',int2str(m),' => vars[i][',int2str(7+k),']'])
%             disp(['    store.impose( new PeqQ(vars[i][',int2str(7+k),'],',variableTemp{1},'[i-1][',int2str((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + (nEtat+nEntre) -1),']) );'])     
        disp(['    store.impose(new IfThen(new PeqC(vars[i-1][5],',int2str(m),') , new PeqQ(vars[i][',int2str(7+k),'],',variableTemp{1},'[i-1][',int2str((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + (nEtat+nEntre) -1),']) ) );'])     
    end
    disp('}')

    disp(['for(int i = 1 ; i < step ; i++){'])
    for k = 1:nSortie % Calcul de Yk = C*Xk + D*Uk

        disp('    // calcul des C(k,:)*Xk')
        for j = 1:nEtat % pour toutes les colonnes de C             
            %disp(((k-1)*(nEtat+nEntre)*2 + j) + nTemp1)
            ind = ((k-1)*(nEtat+nEntre)*2 + j) + nTemp1 -1;
            disp(['    store.impose( new PmulCeqR(vars[i][',int2str(7+j),'],',num2str(C(k,j)),',',variableTemp{1},'[i][',int2str(ind),']) );'])
        end

        disp('    // calcul des D*Uk')
        for j = 1:nEntre             
            %disp((k*nEtat + (k-1)*(nEtat+nSortie+1) + j) + nTemp1) 
            ind = (k*nEtat + (k-1)*(nEtat+nSortie+1) + j) + nTemp1 -1;
            disp(['    store.impose( new PmulCeqR(vars[i][',int2str(j-1),'],',num2str(D(k,j)),',',variableTemp{1},'[i][',int2str(ind),']) );'])            
        end

        disp('    // calcul des sommes (c11*x1k + c12*x2k ...)')
        ind = ((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + 1) + nTemp1 -1;
        disp(['    store.impose( new PeqQ(',variableTemp{1},'[i][',int2str(ind),'],',variableTemp{1},'[i][',int2str(((k-1)*(nEtat+nEntre)*2 + 1) + nTemp1 -1),']) );'])            
        for j = 2:(nEtat+nEntre)            
            %disp(((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + j) + nTemp1)
            ind = ((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + j) + nTemp1 -1;
            disp(['    store.impose( new PplusQeqR(',variableTemp{1},'[i][',int2str(ind-1),'],',variableTemp{1},'[i][',int2str(((k-1)*(nEtat+nEntre)*2 + j) + nTemp1 -1),'],',variableTemp{1},'[i][',int2str(ind),']) );'])            
        end

%             disp(['    // vars[i][',int2str(1+k),'] = C(1,:)*x(k) + D(1,:)*u(k)'])
        disp(['    // m(i) = ',int2str(m),' => vars[i][',int2str(1+k),']'])
        disp(['    store.impose(new IfThen(new PeqC(vars[i][5],',int2str(m),') , new PeqQ(vars[i][',int2str(1+k),'],',variableTemp{1},'[i][',int2str(((k)*(nEtat+nEntre) + (k-1)*(nEtat+nEntre) + (nEtat+nEntre) ) + nTemp1 -1),']) ) );'])            

    end
    disp('}')
    
    
%% Génération des ellipses en contraintes CSP 
clc

x={'Va','Va','alpha'};
y={'h','alpha','h'};
variableTemp = {'temp48Hold','temp45Hold','temp58Hold'};

for k =1:length(variableTemp)
    s = svd(P0{k}'*P0{k});
    disp(['//contrainte pour : (',x{k},',',y{k},') dans le mode HOLD'])

    disp(['FloatVar[] ',variableTemp{k},' = new FloatVar[5];'])
    disp(['for(int i = 0 ; i < 5 ; i++){'])
        disp(['    ',variableTemp{k},'[i] = new FloatVar(store, "',variableTemp{k},'"+i, -1e150, 1e150);']);
    disp('}')

    disp(['store.impose( new PmulQeqR(',x{k},',',x{k},',',variableTemp{k},'[',int2str(0),']) ); // x*x = ',variableTemp{k},'0'])
    disp(['store.impose( new PmulQeqR(',y{k},',',y{k},',',variableTemp{k},'[',int2str(3),']) ); // y*y = ',variableTemp{k},'2'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(0),'],',num2str(s(1)^2),',',variableTemp{k},'[',int2str(1),']) ); // x^2/s1^2 = ',variableTemp{k},'1'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(2),'],',num2str(s(2)^2),',',variableTemp{k},'[',int2str(3),']) ); // y^2/s2^2 = ',variableTemp{k},'3'])
    disp(['store.impose( new PplusQeqR(',variableTemp{k},'[',int2str(1),'],',variableTemp{k},'[',int2str(3),']',',',variableTemp{k},'[',int2str(4),']) ); // x^2/s1^2 + y^2/s2^2 = ',variableTemp{k},'4'])
    disp(['store.impose( new PltC(',variableTemp{k},'[',int2str(4),'],',num2str(1),') ); // x^2/s1^2 + y^2/s2^2 < 1'])
    disp(' ')
    disp(' ')
end

variableTemp = {'temp48Climb','temp45Climb','temp58Climb'};

for k =1:length(variableTemp)
    s = svd(P1{k}'*P1{k});
    disp(['//contrainte pour : (',x{k},',',y{k},') dans le mode CLIMB'])

    disp(['FloatVar[] ',variableTemp{k},' = new FloatVar[5];'])
    disp(['for(int i = 0 ; i < 5 ; i++){'])
        disp(['    ',variableTemp{k},'[i] = new FloatVar(store, "',variableTemp{k},'"+i, -1e150, 1e150);']);
    disp('}')

    disp(['store.impose( new PmulQeqR(',x{k},',',x{k},',',variableTemp{k},'[',int2str(0),']) ); // x*x = ',variableTemp{k},'0'])
    disp(['store.impose( new PmulQeqR(',y{k},',',y{k},',',variableTemp{k},'[',int2str(3),']) ); // y*y = ',variableTemp{k},'2'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(0),'],',num2str(s(1)^2),',',variableTemp{k},'[',int2str(1),']) ); // x^2/s1^2 = ',variableTemp{k},'1'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(2),'],',num2str(s(2)^2),',',variableTemp{k},'[',int2str(3),']) ); // y^2/s2^2 = ',variableTemp{k},'3'])
    disp(['store.impose( new PplusQeqR(',variableTemp{k},'[',int2str(1),'],',variableTemp{k},'[',int2str(3),']',',',variableTemp{k},'[',int2str(4),']) ); // x^2/s1^2 + y^2/s2^2 = ',variableTemp{k},'4'])
    disp(['store.impose( new PltC(',variableTemp{k},'[',int2str(4),'],',num2str(1),') ); // x^2/s1^2 + y^2/s2^2 < 1'])
    disp(' ')
    disp(' ')
end

variableTemp = {'temp48Descent','temp45Descent','temp58Descent'};

for k =1:length(variableTemp)
    s = svd(P2{k}'*P2{k});
    disp(['//contrainte pour : (',x{k},',',y{k},') dans le mode DESCENT'])

    disp(['FloatVar[] ',variableTemp{k},' = new FloatVar[5];'])
    disp(['for(int i = 0 ; i < 5 ; i++){'])
        disp(['    ',variableTemp{k},'[i] = new FloatVar(store, "',variableTemp{k},'"+i, -1e150, 1e150);']);
    disp('}')

    disp(['store.impose( new PmulQeqR(',x{k},',',x{k},',',variableTemp{k},'[',int2str(0),']) ); // x*x = ',variableTemp{k},'0'])
    disp(['store.impose( new PmulQeqR(',y{k},',',y{k},',',variableTemp{k},'[',int2str(3),']) ); // y*y = ',variableTemp{k},'2'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(0),'],',num2str(s(1)^2),',',variableTemp{k},'[',int2str(1),']) ); // x^2/s1^2 = ',variableTemp{k},'1'])
    disp(['store.impose( new PdivCeqR(',variableTemp{k},'[',int2str(2),'],',num2str(s(2)^2),',',variableTemp{k},'[',int2str(3),']) ); // y^2/s2^2 = ',variableTemp{k},'3'])
    disp(['store.impose( new PplusQeqR(',variableTemp{k},'[',int2str(1),'],',variableTemp{k},'[',int2str(3),']',',',variableTemp{k},'[',int2str(4),']) ); // x^2/s1^2 + y^2/s2^2 = ',variableTemp{k},'4'])
    disp(['store.impose( new PltC(',variableTemp{k},'[',int2str(4),'],',num2str(1),') ); // x^2/s1^2 + y^2/s2^2 < 1'])
    disp(' ')
    disp(' ')
end
    
%% Domaines de stabilité global : affichage des ellipses pertinentes

% PROJECTION DANS (Va,h)

figure

subplot(311)
plot(x0(1,:),y0(1,:))
title('Domaines de stabilités du système global : projection dans le plan (Va,h)')
grid on
xlabel('Va')
ylabel('h')
legend('HOLD')

subplot(312)
plot(x1(1,:),y1(1,:),'r')
grid on
xlabel('Va')
ylabel('h')
legend('CLIMB')

subplot(313)
plot(x2(1,:),y2(1,:),'g')
grid on
xlabel('Va')
ylabel('h')
legend('DESCENT')

% PROJECTION DANS (Va,alpha)
figure

subplot(311)
plot(x0(2,:),y0(2,:)*0.005)
title('Domaines de stabilités du système global : projection dans le plan (Va,alpha)')
grid on
xlabel('Va')
ylabel('alpha')
legend('HOLD')

subplot(312)
plot(x1(2,:),y1(2,:)*0.005,'r')
grid on
xlabel('Va')
ylabel('alpha')
legend('CLIMB')

subplot(313)
plot(x2(2,:),y2(2,:)*0.005,'g')
grid on
xlabel('Va')
ylabel('alpha')
legend('DESCENT')

% PROJECTION DANS (alpha,h)
figure

subplot(311)
plot(x0(3,:)*0.005,y0(3,:))
title('Domaines de stabilités du système global : projection dans le plan (alpha,h)')
grid on
xlabel('alpha')
ylabel('h')
legend('HOLD')

subplot(312)
plot(x1(3,:)*0.005,y1(3,:),'r')
grid on
xlabel('alpha')
ylabel('h')
legend('CLIMB')

subplot(313)
plot(x2(3,:)*0.005,y2(3,:),'g')
grid on
xlabel('alpha')
ylabel('h')
legend('DESCENT')

%% Gestion de l'automate

%Declaration des etats et initialisation
u;

EP = 0; % Etat present
ES = EP; % Etat suivant

sys_hold = sys_holdk5;
sys_climb = sys_climbk5;
sys_descent = sys_descentk5;


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
