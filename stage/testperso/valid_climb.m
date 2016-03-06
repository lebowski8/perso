%clear all
home;
% load('../modele.mat')
load('modele.mat')

% X(:,1) = [0 0 0 235 0 0 0 1000 0 0 0]';
  
% etat stable pour hold
X(:,1) = 1e3*[0.0003
    0.0001
   -0.0000
    0.2358
   -0.0282
   -0.0000
   -0.0001
    1.0000
    0.1147
    0.0001
    0.0000]';


Ts=0.5;                % Pas de temps en seconde
tsim=300*Ts;                % Temps de simulation
t = 0:Ts:tsim;       % Vecteur de temps

u=[235*ones(1,size(t,2)); 1000*ones(1,size(t,2))];
cpt=1;
sys = sys_hold;

for k=1:numel(t)
    
    X(:,cpt+1)=sys.A*X(:,cpt)+sys.B*u(:,cpt);
    Y(:,cpt)=sys.C*X(:,cpt)+sys.D*u(:,cpt);
    
%         disp('MODE DESCENT')
%         disp('stable à Xk ?')
        if(inEllipse(P_hold,X(:,cpt)))
            disp('stable')
        end
%         disp('stable à Xk+1 ?')
%         inEllipse(P_hold,X(:,cpt+1))    
    
    Yk=Y(:,cpt);
    h(cpt) = Yk(2);
    cpt=cpt+1;
end

figure;plot(t,h)
%figure;step(sys_climb)