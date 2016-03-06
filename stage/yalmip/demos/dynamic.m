


%% premier exemple : http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Basics
clear all
close all
clc

% Define variables
x = sdpvar(2,1);
% Define constraints and objective
Constraints = [sum(x) <= 1, x(1)==0, x(2) >= 0.5];
Objective = x'*x+norm(x);
% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','sedumi');
% Solve the problem
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution = value(x)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

%% deuxième exemple : Linear systems and Lyapunov stability (2) commande_etat1.pdf
clear all
close all
clc

% définition des variables
A = [0 1 ; 1 5];
P = sdpvar(2,2, 'symmetric');

% définition des contraintes
F = P>0;
F = F + ((-[eye(2) A']*[zeros(2) P ; P zeros(2)]*[eye(2) ; A]) > 0)

% définition des options du solver
options = sdpsettings('verbose',1,'solver','sedumi');

% résolution du problème par le solver
sol = solvesdp(F,trace(P),options)

% affichage de P
value(P)

[vecteur valeur] = eig(value(P))

 x = [-5:1:5 ; -5:1:5];
 
 for i = 1:length(x)
     V(i) = x(:,i)'*value(P)*x(:,i);
     %dV = x'*A'*value(P)*x + x'*value(P)*A*x;
 end
 
xx1=-5:1:5;
xx2=-5:1:5;
[XXX,YYY]=meshgrid(xx1,xx2);
ZZZ=XXX.^2*1+YYY.^2*3;
figure
meshc(XXX,YYY,ZZZ)
figure
[H,c] = contour(XXX,YYY,ZZZ);

% plan de coupe : hslice = slice(D,50,[],[]);


%% exemple 2 bis
clear all
close all
clc

% définition des variables
A = [0 1 ; -1 -2];
P = sdpvar(2,2, 'symmetric');
x = sdpvar(2,1);

% définition des contraintes
F = P>0;
F = F + ((-[eye(2) A']*[zeros(2) P ; P zeros(2)]*[eye(2) ; A]) > 0)

% définition des objectifs
Objective = eig(P);

% définition des options du solver
options = sdpsettings('verbose',1,'solver','sedumi');

% résolution du problème par le solver
sol = solvesdp(F,trace(P),options)

% affichage de P
value(P)

vp = eig(value(P))

 V = x'*value(P)*x;
% dV = dx'*P*x + x'*P*dx
 dV = x'*A'*value(P)*x + x'*value(P)*A*x;

% test d'affichage
plot(F,[x(1),x(2)],[],[],options)

%% exemple 3
clear all
close all
clc

x = sdpvar(3,1);
F = [
    [1 x(2);x(2) x(1)] >= 0,
    3 >= x,
    [1 x';x eye(3)] >= 0,
    sum(x) >= x'*x
    ]

plot(F)


