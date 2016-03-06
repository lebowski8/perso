clear all
close all
home





%% Modèle Longitudinal : représentation d'état xdot = Ax + Bu, y = Cx
% y = [V alpha theta q gamma]'
A =[-0.018223 -0.088571 -9.78 0
    -0.003038 -1.2563 0  1
    0 0 0 1
    0.0617 -28.075 0 -4.593];
B  = [ 0 1.1962
    0 -0.0012 
    0 0
    -7.84 -4.05]; % /!\ Changement par rapport à l'énoncé (signe opposé pour modèle cohérent)
C= [1 0 0 0
    0 57.296 0 0
    0 0 57.296  0
    0 0 0 57.296
    0 -57.296 57.296 0];
D=zeros(5,2);
Slong = ss(A,B,C,D);


%% Stabilité au sens de Lyapunov
P = sdpvar(4,4, 'symmetric');

% définition des contraintes
F = P>0;
%F = F + ((-[eye(2) A']*[zeros(2) P ; P zeros(2)]*[eye(2) ; A]) > 0);
F1 = A'*P + P*A;
F = F + (F1<0);
% définition des options du solver
options = sdpsettings('verbose',1,'solver','sedumi');

% résolution du problème par le solver
sol = solvesdp(F,[],options);
checkset(F);

% on vérifie que P est définie positive
valP = value(P);
spectreP = eig(valP);
flagP = 0;
for i = 1:rank(valP)
    if spectreP(i) <= 0 
        disp(spectreP(i))
        flagP = 1;
    end
end

if flagP == 1
    disp('P : non définie positive')
else
    disp('P : définie positive')
end



%% Calcul de V(x) et dV(x)s 

x1=-10:1:10;
x2=-10:1:10;
x3=-10:1:10;
x4=-10:1:10;

for cpt1=1:numel(x1)
    for cpt2=1:numel(x2)
        for cpt3=1:numel(x3)
            for cpt4=1:numel(x4)
                x = [x1(cpt1) x2(cpt2) x3(cpt3) x4(cpt4)];
                myV=x*valP*x';
                myLyap(cpt1,cpt2,cpt3,cpt4) = myV; 
                if myV<0
                    ind_error=[cpt1 cpt2 cpt3 cpt4] 
                end
            end
        end
    end
end

for cpt1=1:numel(x1)
    for cpt2=1:numel(x2)
        for cpt3=1:numel(x3)
            for cpt4=1:numel(x4)
                x = [x1(cpt1) x2(cpt2) x3(cpt3) x4(cpt4)];
                mydV=x*A'*valP*x' + x*valP*A*x';
                mydLyap(cpt1,cpt2,cpt3,cpt4) = mydV; 
            end
        end
    end
end


%% Vérification
if(all(all(myLyap>=0)))
    disp('V(x) est positive pour tout x')
end

if(all(all(mydLyap<=0)))
    disp('dV(x)/dx est négative pour tout x')
end

%% projection de V(x) sur les différents plans
p12 = squeeze(myLyap(:,:,1,1));
p13 = squeeze(myLyap(:,1,:,1));
p14 = squeeze(myLyap(:,1,1,:));
p23 = squeeze(myLyap(1,:,:,1));
p24 = squeeze(myLyap(1,:,1,:));
p34 = squeeze(myLyap(1,1,:,:));


%% Affichage 

figure
subplot(211)
h = meshc(x1,x2,p12);
title('V(x)')
xlabel('X(1)')
ylabel('X(2)')
zlabel('V(X)')
subplot(212)
contour(x1,x2,p12);
xlabel('X(1)')
ylabel('X(2)')

figure
subplot(211)
h = meshc(x1,x3,p13);
title('V(x)')
xlabel('X(1)')
ylabel('X(3)')
zlabel('V(X)')
subplot(212)
contour(x1,x3,p13);
xlabel('X(1)')
ylabel('X(3)')

figure
subplot(211)
h = meshc(x1,x4,p14);
title('V(x)')
xlabel('X(1)')
ylabel('X(4)')
zlabel('V(X)')
subplot(212)
contour(x1,x4,p14);
xlabel('X(1)')
ylabel('X(4)')

figure
subplot(211)
h = meshc(x2,x3,p23);
title('V(x)')
xlabel('X(2)')
ylabel('X(3)')
zlabel('V(X)')
subplot(212)
contour(x2,x3,p23);
xlabel('X(2)')
ylabel('X(3)')

figure
subplot(211)
h = meshc(x2,x4,p24);
title('V(x)')
xlabel('X(2)')
ylabel('X(4)')
zlabel('V(X)')
subplot(212)
contour(x2,x4,p24);
xlabel('X(2)')
ylabel('X(4)')

figure
subplot(211)
h = meshc(x3,x4,p34);
title('V(x)')
xlabel('X(3)')
ylabel('X(4)')
zlabel('V(X)')
subplot(212)
contour(x3,x4,p34);
xlabel('X(3)')
ylabel('X(4)')


% figure(1)
% [C h] = contour(x1,x2,myLyap);
% title('projection de V(x)')
% hold on
% grid on
% xlabel('X(1)')
% ylabel('X(2)')

% hold on
% 
% plot(C(1,2:34),C(2,2:34),'r')
% plot(C(1,36:84),C(2,36:84),'g')
% plot(C(1,86:146),C(2,86:146),'black')




