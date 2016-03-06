clear all
%close all
%clc
home

% Cf GOMES & Tarbouriech 2005 
% POUR MODEL CONTINUE

% n = 2;
% c = 2;

% 
% xn = zeros(n,1);
% xc = zeros(c,1);

%% exemple 1
% le système
sys = ss(); 

sys.a = 0.1;
sys.b = 1;
sys.c = 1;
sys.d = 0;

n = length(sys.a)

% le contrôleur
Ac = 0;
Bc = -0.2;
Cc = 1;
Dc = -2;

sysc = ss(Ac,Bc,Cc,Dc);

dimAc = size(sysc.a);
c = length(sysc.a)

% définition de la zone
uo = 1; % control bound 
m = length(uo)
Theta0 = {[1;1] , [1;-1] , [-1;1] , [-1;-1]}; % défini la zone à respecter
r = length(Theta0);

%% exemple 2
% le système
% sys = ss(); 
% 
% sys.a = [0.1 -0.1 ; 0.1 -3];
% sys.b = [5 0 ; 0 1];
% sys.c = [1 0 ; 0 1];
% sys.d = 0;
% 
% n = length(sys.a)
% 
% % le contrôleur
% Ac = [-171.2 27.2 ; -68 -626.8];
% Bc = [-598.2 5.539 ; -4.567 149.8];
% Cc = [0.146 0.088 ; -6.821 -5.67];
% Dc = [0 0 ; 0 0];
% 
% sysc = ss(Ac,Bc,Cc,Dc);
% 
% c = length(Ac)
% 
% % définition de la zone
% uo = [5 ; 2]; % control bound 
% m = length(uo)
% Theta0 = {[1;1;0;0] , [1;-1;0;0] , [-1;1;0;0] , [-1;-1;0;0]}; % défini la zone à respecter
% r = length(Theta0);


%% le système bouclé

A = [sys.a + sys.b*Dc*sys.c , sys.b*Cc ; Bc*sys.c , Ac];
B = [sys.b ; zeros(c,m)]; 
R = [zeros(n,n) ; ones(c,n)];
K = [Dc*sys.c , Cc];



%% problème de minimisation (25)
% les dimensions sont détaillées "Corollary 1"
W = sdpvar(n+c,n+c,'symmetric');
Z = sdpvar(n,m);
S = sdpvar(m,m,'diagonal');
Y = sdpvar(m,n+c);
mu = sdpvar(1);



%% définition des contraintes
cdt1 = [W*A'+A*W B*S+R*Z-Y' ; S*B'+Z'*R'-Y -2*S];

F = (cdt1 < 0); 

for i=1:m
    F = F + ([W W*K(i,:)'-Y(i,:)' ; K(i,:)*W-Y(i,:) uo(i)^2] >= 0);
end

for i=1:r
   F = F + ([mu Theta0{i}' ; Theta0{i} W] >= 0); 
end


% définition des options du solver
options = sdpsettings('verbose',1,'solver','sedumi');

% résolution du problème par le solver
sol = solvesdp(F,[mu],options);
checkset(F);

beta = 1/sqrt(value(mu))
Ec = value(Z*inv(value(S)))
P = inv(value(W))


%% ellipse
 close all

% valeur de Tarbouriech
% beta = 1.9165
% P = [0.0497 -0.0377 ; -0.0377 0.1472]
% Ec = 0.092


[v,d] = eig(P);
spectreP = diag(d);
%[x,y] = calculateEllipse(0,0,max(spectreP),min(spectreP),subspace(v(:,2),[1 0]'));


figure
[x y] = projellisa(P,[],'r');

% figure
% plotellisa(P,zeros(length(P),1))
% grid on
% axis equal


hold on
plot(x,y,'r')
axis equal
grid on

% for i=1:length(Theta0)
%     p = beta*Theta0{i};
%     plot(p(1),p(2),'bo')
% end


