clear all
close all
home




A = [0 1 ; -1 -2];
B = [0 ; 0];
C = [1 0 ; 0 1];
D = [0 ; 0];
sys0c=ss(A,B,C,D);

P = sdpvar(2,2, 'symmetric');

% définition des contraintes
F = P>0;
F1 = A'*P + P*A;
F = F + (F1<0);



% définition des options du solver
options = sdpsettings('verbose',1,'solver','lmilab');

% résolution du problème par le solver
sol = solvesdp(F,[trace(P)],options);
checkset(F);

%% on vérifie que P est définie positive (donc V positive pour tout x)
valP = value(P)
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

%% on vérifie que dV est négative pour tout x 
x1=-10:1:10;
x2=-10:1:10;

for cpt1=1:numel(x1)
    for cpt2=1:numel(x2)
   myV=[x1(cpt1) x2(cpt2)]*valP*[x1(cpt1);x2(cpt2)];
   myLyap(cpt1,cpt2) = myV; 
   if myV<0
      ind_error=[cpt1 cpt2] ;
   end
    end
end

for cpt1=1:numel(x1)
    for cpt2=1:numel(x2)
   mydV=[x1(cpt1) x2(cpt2)]*A'*valP*[x1(cpt1);x2(cpt2)] + [x1(cpt1) x2(cpt2)]*valP*A*[x1(cpt1);x2(cpt2)];
   mydLyap(cpt1,cpt2) = mydV; 
    end
end


if(all(all(myLyap>=0)))
    disp('V(x) est positive pour tout x')
end

if(all(all(mydLyap<=0)))
    disp('dV(x)/dx est négative pour tout x')
end

%% calcul du domaine de stabilité (ellipse)
valP = value(P);
[x y] = projellisa(value(P),[],'r');

figure
plot(x, y)
grid on
axis square
xlabel('X(1)')
ylabel('X(2)')
title('Domaine de stabilité, xT*P*x < 1')