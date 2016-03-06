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
tsim=8*Ts;                % Temps de simulation
t = 0:Ts:tsim-Ts;       % Vecteur de temps
i=1;

eps = 50; % variation d'altitude autour du mode maintien
h(1) = 10000*0; % altitude initiale
Va = 235*0; % vitesse initiale
X(:,1) = zeros(nEtat,1);
X(4,1) = Va;
X(8,1) = h(1);




u = [ 235*ones(1,length(t)); 60*ones(1,length(t))];
% on est pas dans l'ellipse formée par le plan (Va,alpha)
% u = [ 235*ones(1,length(t)); linspace(200,300,length(t)/2) 300*ones(1,length(t)/2)];
% u = [ 20*ones(1,length(t)); 200*ones(1,length(t)/2) 100*ones(1,length(t)/2)];

% u = [235*ones(1,length(t)) ; 1499.9 1499.9 1499.9 1126.9 1499.9]; % test planif

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
    
    %% Définition des systèmes et étude de stabilité : 

    P = calcMatriceLyap(sys_hold, sys_descent, sys_climb); % avec les idées de emmanuel

    
%% Etude de la stabilité des système k5

% HOLD EP = 0 :
    [P0 x0 y0 P_holdk5] = domaineLyap(sys_holdk5,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD
%     [P0k10 x0k10 y0k10] = domaineLyap(sys_holdk10,'discret'); % P0 = {P_48 ; P_45 ; P_58} pour le mode HOLD

% CLIMB EP = 1 : 
    [P1 x1 y1 P_climbk5] = domaineLyap(sys_climbk5,'discret');
% DESCENT EP = -1 : 
    [P2 x2 y2 P_descentk5] = domaineLyap(sys_descentk5,'discret');
    
%% Génération des Systèmes pour le format de donnée de Cplex sous JAVA
clc    
    
    mode = {-1, 0, 1};
    sysName = {'_Descent', '_Hold', '_Climb'};
    modes = {sys_descentk5, sys_holdk5, sys_climbk5};
%     modes = {sys_descent, sys_hold, sys_climb}; 
    
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
            disp(['static double[][] ',matName{m}, sysName{k}, ' = {'])
            M = mat{m};
            for i = 1:size(M,1)
                chaine = '{';
                for j = 1:(size(M,2)-1)
                    chaine = [chaine , num2str(M(i,j)),' , '];
                end
                disp(['         ',chaine, num2str(M(i,size(M,2))), '},'])
            end
            disp('};')
        end

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
% plot(x0(2,:),y0(2,:)*0.005)
plot(x0(2,:),y0(2,:))
title('Domaines de stabilités du système global : projection dans le plan (Va,alpha)')
grid on
xlabel('Va')
ylabel('alpha')
legend('HOLD')

subplot(312)
% plot(x1(2,:),y1(2,:)*0.005,'r')
plot(x1(2,:),y1(2,:),'r')
grid on
xlabel('Va')
ylabel('alpha')
legend('CLIMB')

subplot(313)
% plot(x2(2,:),y2(2,:)*0.005,'g')
plot(x2(2,:),y2(2,:),'g')
grid on
xlabel('Va')
ylabel('alpha')
legend('DESCENT')

% PROJECTION DANS (alpha,h)
figure

subplot(311)
% plot(x0(3,:)*0.005,y0(3,:))
plot(x0(3,:),y0(3,:))
title('Domaines de stabilités du système global : projection dans le plan (alpha,h)')
grid on
xlabel('alpha')
ylabel('h')
legend('HOLD')

subplot(312)
% plot(x1(3,:)*0.005,y1(3,:),'r')
plot(x1(3,:),y1(3,:),'r')
grid on
xlabel('alpha')
ylabel('h')
legend('CLIMB')

subplot(313)
% plot(x2(3,:)*0.005,y2(3,:),'g')
plot(x2(3,:),y2(3,:),'g')
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

%% Définition de la mission
coef = 200;
eps = 50;

Hgoal = 4*coef;
Lgoal = 5875;

Hmax = Hgoal+coef;
Lmax = Lgoal+1000;

% zone1 = [2 2 6 6 ; 5 8 8 5]'*coef;
% zone2 = [6 6 6.5 6.5 ; 0 2 2 0]'*coef;
% zone3 = [7 7 8.5 8.5 8 8; 3 4 4 2 2 3]'*coef;
% zone3 = [7 7 8.5 8.5; 3 4 4 3]'*coef;
% zone4 = [9.5 9.5 11.5 11.5 ; 0 2 2 0]'*coef;

zone1 = [2 2 10 10 ; 5 8 8 5]'*100;
zone2 = [30 30 50 50 ; 0 2 2 0]'*100;
% zone3 = [7 7 8.5 8.5 8 8; 3 4 4 2 2 3]'*coef;
zone3 = [20 20 30 30; 4 6 6 4]'*100;
% zone4 = [9.5 9.5 11.5 11.5 ; 0 2 2 0]'*100;

figure
hold on

fill(zone1(:,1),zone1(:,2),'r')
% fill(zone2(:,1),zone2(:,2),'b')
% fill(zone3(:,1),zone3(:,2),'g')
% fill(zone4(:,1),zone4(:,2),'y')

axis([-coef Lmax+coef -coef Hmax+coef])
% axis equal
xlabel('distance (km)')
ylabel('hauteur (km)')
grid on

text(0,-0.5*coef,[' Départ']);
plot(0, 0,'r.','MarkerSize',20) % On affiche les points
text(Lgoal+0.5*coef,Hgoal,['  Arrivé']);
plot(Lgoal, Hgoal,'r.','MarkerSize',20) % On affiche les points

rectangle('Position',[Lgoal-eps, Hgoal-eps, eps*2, eps*2]) % on affiche la zone de tolérance à l'arrivée
text(Lgoal-eps-350,Hgoal-eps+150,['  Zone de tolérance']);

%% Affichage des résultats Cplex Java
% t = [0.0 2.5 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 27.5 ];
% Va = [0.0 35.493567955376946 78.66924390243786 134.69817423267008 194.1775442128696 248.0106195959816 290.0864431333358 317.8273457781422 335.71018258157096 345.7736720492155 349.55320655839614 350.0 ];
% h = [0.0 4.802707141699997 23.28332505028319 126.50413537607508 263.85897649123484 376.4720108573671 454.5333036627975 547.4206393278015 629.1161451254357 688.6883111558681 729.0727128417208 750.0 ];
% alpha = [2.7000623958883807E-12 -0.016536607591216532 -0.044470184340699326 -0.04984088924609312 -0.1981268223197943 -0.24798870702694958 -0.27603583271971355 -0.2196725528696591 -0.3251604467485786 -0.32532386581989386 -0.3221731677146714 -0.32148054759964484 ];
% L = [0.0 0.0 88.73391988844236 285.40702964453703 622.1524652262123 1107.5963257583862 1727.62287474834 2452.8389825816794 3247.4073470270355 4086.682803480963 4951.116983604002 5824.999999999992 ];
% x0 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]';
% x1 = [18.837048988711704 0.002259577890056672 0.002696769283032035 8.460786466061478 -0.29518266190396897 -1.3620892989099856E-5 1.255101319127694E-4 0.33170646865440445 37.04670231879328 9.421070549485162E-6 0.0 ]';
% x2 = [35.82181610540525 0.011449205596424351 0.0066945810460367055 41.47049031447705 -1.9113213408326932 0.0032020593079790257 0.0050703923287351476 5.3329605372716635 60.32507021432548 -0.005761275493350688 -0.6366488152216717 ]';
% x3 = [49.31866724146468 -0.036306930324098326 -0.0012532031956625635 87.9963932258327 12.940929671092407 0.09499935531007181 0.13329327042577568 29.180925050283186 95.98051970314845 -0.15797467676193364 -12.538546270214548 ]';
% x4 = [56.880676164049646 0.06617353694217366 0.030984620922963207 144.83614132963936 -8.356339845860543 -0.05711605108829014 0.18419812563830956 126.86004780322348 125.15187759751153 -0.09617062045429337 0.0 ]';
% x5 = [54.374238149555026 0.12677625333257736 0.018550563490491834 205.011128786393 -35.81337663754232 -0.04438424394980456 0.046269894812743795 265.5173007851877 147.05157993846757 -0.015359475117349947 -0.42704465407428316 ]';
% x6 = [44.767434645068704 0.15480316283925172 0.009035149436405604 259.05419482751745 -43.81993075648347 -0.037115288808832914 -0.05384221999647143 378.6529542258803 161.80622634155145 0.05448765309533174 -0.45080484193174186 ]';
% x7 = [33.44388170880852 0.12302103410706854 -0.0033122559686777527 300.75143547696126 -34.702516619649444 0.03594617013256429 -0.052900099437732745 460.4309036627974 170.55478341207777 0.007459038677521568 -7.439998621337826 ]';
% x8 = [22.860036007664316 0.1697254481347582 0.011905568894942601 328.94267275043586 -42.71682885039276 -0.042155270484708116 -0.05847049736005312 547.7337829459284 174.95384438115906 0.06413401275464889 0.0 ]';
% x9 = [14.241793973139693 0.1826452208049381 0.0031168498423924723 346.9136526312998 -51.674541454471246 -0.018342896284154718 -0.1391527331673593 633.2430383192669 176.5500002596642 0.10302754531936571 -1.0490301581720096 ]';
% x10 = [7.437783246596155 0.1811780986713037 -6.492966994749994E-4 356.82876968753453 -50.16365802978174 -0.011330076062391559 -0.17398617171651942 693.1401481669707 175.99602400903473 0.1264613906510492 -1.9394498000712224 ]';
% x11 = [2.717305849655251 0.18387337523685904 6.8678891014431855E-6 360.73377966098974 -50.622066252299824 -0.015586428415324267 -0.20803477732451686 733.1993774039242 175.03829137766638 0.15447406763821192 -1.511115543054443 ]';

t = [0.0 2.5 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 27.5 ];
Va = [0.0 41.43338906104834 104.35582724317362 169.61491747855618 206.05734055926996 242.05154034444814 274.5857609662401 301.14198849624063 319.2798360326799 331.9359876482385 339.5434121701044 324.120571035051 ];
h = [0.0 5.592264459706881 24.289417817330722 129.2490282519766 267.15565654601346 398.9134795477539 494.83690449501563 552.9488706940924 620.1551699170228 677.1488827116484 727.4746031573493 750.0 ];
alpha = [2.7000623958883807E-12 -0.01931112170553838 -0.06080070822356375 -0.07121980761685318 -0.2126587139558725 -0.20807430909854366 -0.2819179380430896 -0.28925357958997194 -0.21680254184501765 -0.31542093373851077 -0.24834292784218676 -0.3114290416283985 ];
L = [0.0 0.0 103.58347265262084 364.4730407605549 788.5103344569454 1303.6536858551203 1908.7825367162407 2595.246939131842 3348.1019103724434 4146.301500454143 4976.141469574739 5825.0 ];
x0 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]';
x1 = [21.988879112417003 0.0026521217249 0.00315017725546 9.876751918962 -0.34883200847 -3.618442439999997E-5 1.195877819E-4 0.38609836673 43.24600397567049 4.311233839999999E-5 0.0 ]';
x2 = [48.290054097822214 0.018874052202182235 0.009596891057979775 51.413166831776145 -3.684531110519218 -0.0028228733109457883 -0.002580390000403894 5.980141046187231 83.34838836293345 0.0036645470086336058 0.009978325879313427 ]';
x3 = [61.28932336686525 -0.024855188607612027 9.466035251041779E-4 113.77357636359793 9.866428959246065 0.09488393240249507 0.11781342352275015 30.187017817330727 116.78450148209734 -0.14716600151287695 -11.881590574444536 ]';
x4 = [48.44565605278924 0.08030174379780695 0.030597416294249353 174.2295168409048 -12.08389251935366 -0.05763564679696254 0.16841254019642535 129.40126439121227 117.39900589234061 -0.08478402346671393 0.0 ]';
x5 = [37.67649740235677 0.11214128486750667 0.010955391787002718 216.72522860425318 -32.08166384744125 -0.0118558185865817 0.07133786310405914 270.5520555555694 136.63367825866953 -0.054872218929087016 -4.186366398115093 ]';
x6 = [33.61292050924579 0.15299711175968164 0.008947272276815499 253.11956040662483 -41.9815285878746 -0.049898277049377146 -0.005968277018047051 399.24641471506453 151.07579039269885 0.027257243440617966 0.0 ]';
x7 = [28.004973131564583 0.16750525223275115 0.0039946056652046865 285.73289419824107 -48.37756207172521 -0.03150500710976522 -0.10903227705837593 497.67100901960004 161.23763424864802 0.09004737987803652 0.08155432572638323 ]';
x8 = [21.8210907319894 0.13252603723213466 -0.00568522693940464 311.8837200854738 -37.44884935806536 0.03150676054524658 -0.10569168509408514 558.8464706940923 167.70173696759372 0.0493458084513463 -5.590361588895683 ]';
x9 = [15.768043294144622 0.1653515800774663 0.008441343704746924 330.3785973795717 -42.055421161647914 -0.031174684355227564 -0.10309344357798811 620.458620204258 171.37223977621434 0.09019712839593787 0.0 ]';
x10 = [10.700305696863667 0.1699903351824697 0.001447704111976306 343.0062804986943 -47.50490827588091 -0.0067347087740791504 -0.15439232754461188 681.7979134408298 173.14815303547874 0.10813222732084243 -1.4659994702467098 ]';
x11 = [-4.037557917190195 0.17730733225865794 -6.60757546685774E-4 346.11534769960355 -48.64638677822348 -0.01835080066798588 -0.1853623853900228 727.5826701561133 153.0265844464836 0.1395289031130047 0.0 ]';

% figure
% plot(t,Va)
% title('Va(t)')
% grid on
% hold on
% 
% figure
% plot(t,alpha)
% title('alpha(t)')
% grid on
% hold on

hold on
plot(L,h,'b')
title('Déplacement longitudinal')
xlabel('Distance (m)')
ylabel('Altitude (m)')
grid on
hold on


%% Verif stabilité avec énergie décroissante
% 0 : maintien
% 1 : descent
% 2 : monté

P00 = P{1};
P01 = P{2};
P02 = P{3};

P11 = P{4};
P10 = P{5};
P12 = P{6};

P22 = P{7};
P20 = P{8};
P21 = P{9};


V(1) = calcLyap(P20,x0);
V(2) = calcLyap(P00,x1);
V(3) = calcLyap(P01,x2);
V(4) = calcLyap(P10,x3);
V(5) = calcLyap(P00,x4);
V(6) = calcLyap(P01,x5);
V(7) = calcLyap(P10,x6);
V(8) = calcLyap(P00,x7);
V(9) = calcLyap(P01,x8);
V(10) = calcLyap(P11,x9);

figure
plot(t,V)




