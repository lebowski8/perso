clear all
close all
%clc
home

% Cf commande_etat1.pdf page 4/32


% définition des variables
A = [0 1 ; -1 -2];
%A = [0 4 ; -8 -12];
%A = [0 1 0 ; -1 -2 0 ; 0 1 -1];
P = sdpvar(2,2, 'symmetric');
Q = sdpvar(2,2, 'symmetric');

% sys = ss();
% sys.a = A; 
% sys.b=[0;1];
% sys.c = [1 0];

% définition des contraintes
F = P>0;
F1 = A'*P + P*A;
F = F + (F1<0);

% F = (Q > 0) + (Q*A' + A*Q < 0);

% définition des options du solver
options = sdpsettings('verbose',1,'solver','lmilab');

% résolution du problème par le solver
sol = solvesdp(F,[trace(P)],options);
checkset(F);


% on vérifie que P est définie positive
valP = value(P);
% valP = inv(value(Q));
valP
spectreP = eig(valP);
flagP = 0;
for i = 1:rank(valP)
	if spectreP(i) <= 0 
        disp(spectreP(i));
        flagP = 1;
    end
end

if flagP == 1
	disp('P : non définie positive')
else
	disp('P : définie positive')
end

% hold on
% [x y] = projellisa(value(P),[],'r');
% bleu : -trace() rouge : []
% green : trace(Q), bleu : trace(P)
% plot(x, y)

%% test d'affichage comme dans le PDF page 32 de commade_etat1
% home
% clear p11 p12 p22 dV
% syms p11 p12 p22
% %syms p22
% %p11 = 0;
% %p12 = 1;
% dV = [-2*p12 p11-2*p12-p22 ; p11-2*p12-p22 2*p12-4*p22];
% vp = eig(dV);
% s1 = solve(vp(1),p22);
% s2 = solve(vp(2),p22);
% 
% p221 = s1(1);
% p222 = s1(2);
% 
% echant = 0:0.1:1;
% 
% for cpt1 = 1:numel(echant)
%     for cpt2 = 1:numel(echant)
%         p11 = echant(cpt1);
%         p12 = echant(cpt2);
%         subs(p221);
%         test(cpt1,cpt2) = abs(subs(p221));
%     end
% end
% 
% figure
% meshc(echant,echant,test);
% title('p22')
% xlabel('p11')
% ylabel('p12=p21')
% zlabel('p22')

 %%
 %clc
 %close all
 
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


%% Vérification
if(all(all(myLyap>=0)))
    disp('V(x) est positive pour tout x')
end

if(all(all(mydLyap<=0)))
    disp('dV(x)/dx est négative pour tout x')
end
%% Affichage 

[X,Y]=meshgrid(x1,x2);

figure
h = meshc(X,Y,myLyap);
title('V(x)')
xlabel('X(1)')
ylabel('X(2)')
zlabel('V(X)')

% figure
% mesh(X,Y,mydLyap)
% title('V(x)/dx')
% xlabel('X(1)')
% ylabel('X(2)')
% zlabel('V(X)')

% figure
% [C h] = contour(X,Y,myLyap);
% clabel(C,h);
% title('projection de V(x)')
% grid on

[x y] = projellisa(value(P),[],'r');
figure
plot(x, y)
grid on
axis square
xlabel('X(1)')
ylabel('X(2)')


% figure
% C1 = contourc(x1,x2,myLyap,[1 1]);
% plot(C1(1,2:end),C1(2,2:end))
% title('projection de V(x)')
% grid on
% xlabel('X(1)')
% ylabel('X(2)')


% [v,d] = eig(valP)
% [x,y] = calculateEllipse(0,0,spectreP(2),spectreP(1),0);
% hold on
% %plot(x,y,'r')
% plot(v(1,1),v(2,1),'ro')
% plot(v(1,2),v(2,2),'g+')


% plan de coupe : hslice = slice(D,50,[],[]);





