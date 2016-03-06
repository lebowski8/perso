function [P_xy x0 y0 valP] = domaineLyap(sys,mode)
%DOMAINELYAP
% premet de tester la stabilité du système et de sortir le domaine de
% stabilité : x'Px < 1
% mode : 'continu' ou 'discret' (continu par défaut)
%

A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;
%A = [0 1 0 ; -1 -2 0 ; 0 1 -1]
%A = [0 5 ; -1 -2];

P = sdpvar(length(A),length(A), 'symmetric');
% s = sdpvar(1,1);
% s1 = sdpvar(1,1);
% Q = sdpvar(2,2, 'symmetric');

C = zeros(size(C));
D = zeros(size(D));
if nargin > 1
   if mode == 'discret'
        F = P>0;
         F1 = A'*P*A - P;
% Permet la prise en compte de B*uk en discret ! Cf Analysis of switched
% normal discrete-time systems
%         F1 = [(A'*P*A-P+C'*C) , (A'*P*B+C'*D) ; (B'*P*A+D'*C) , (B'*P*B-eye(size(B,2))+D'*D)];
        F1 = [(A'*P*A-P+C'*C) , (A'*P*B+C'*D) ; (B'*P*A+D'*C) , (B'*P*B-eye(size(B,2))+D'*D)];
        F = F + (F1<0);   
        
%         s1 = eig(P);
%         s = max(s1);
   else
        F = P>0;
        F1 = A'*P + P*A;
        F = F + (F1<0);
   end
    
else
    F = P>0;
    F1 = A'*P + P*A;
    F = F + (F1<0);
end
% définition des contraintes


% Q = inv(P) % permet de maximiser P
% F = (Q > 0) + (Q*A' + A*Q < 0);

% définition des options du solver
options = sdpsettings('verbose',1,'solver','lmilab');
% options = sdpsettings('verbose',1,'solver','sedumi');

% résolution du problème par le solver
sol = solvesdp(F,[],options);
% sol = solvesdp(F,[trace(-P)],options); 
% sol = solvesdp(F,[max(eig(P))],options); % j'aurai voulu minimiser la
% valeur max des valeurs propres (singulières) de P, pour avoir de "jolie"
% ellipses...
% sol = solvesdp(F,[(trace(-P'*P))],options);
% sol = solvesdp(F,[trace(Q)],options);
checkset(F);
% P = inv(value(Q))

%% on vérifie que P est définie positive (donc V positive pour tout x)
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

%% on vérifie que dV est négative pour tout x 
% x1=-10:1:10;
% x2=-10:1:10;
% 
% x(:,1) = -10:1:10;
% x(:,2) = -10:1:10;
% 
% for cpt1=1:numel(x1)
%     for cpt2=1:numel(x2)
%    myV=[x1(cpt1) x2(cpt2)]*valP*[x1(cpt1);x2(cpt2)];
%    myLyap(cpt1,cpt2) = myV; 
%    if myV<0
%       ind_error=[cpt1 cpt2] ;
%    end
%     end
% end
% 
% 
% 
% myLyap2 = x*valP*x';
% 
% for cpt1=1:numel(x1)
%     for cpt2=1:numel(x2)
%    mydV=[x1(cpt1) x2(cpt2)]*A'*valP*[x1(cpt1);x2(cpt2)] + [x1(cpt1) x2(cpt2)]*valP*A*[x1(cpt1);x2(cpt2)];
%    mydLyap(cpt1,cpt2) = mydV; 
%     end
% end
% 
% mydLyap2 = x*A'*valP*x' + x*valP*A*x';
% 
% 


% if nargin > 1
%    if mode == 'discret'
%         if(all(eig(A'*valP*A - valP) >= 0))
%             disp('AtPA-P est définie positive')
%         end
% 
%         if(all(eig(A'*valP*A - valP) <= 0))
%             disp('AtPA-P est définie négative')
%         end
%    else
%         if(all(eig(A'*valP + valP*A) >= 0))
%             disp('AtP + PA est définie positive')
%         end
% 
%         if(all(eig(A'*valP + valP*A) <= 0))
%             disp('AtP + PA est définie négative')
%         end
%    end
%     
% else
%     if(all(eig(A'*valP + valP*A) >= 0))
%         disp('AtP + PA est définie positive')
%     end
% 
%     if(all(eig(A'*valP + valP*A) <= 0))
%         disp('AtP + PA est définie négative')
%     end
% end


%% calcul du domaine de stabilité (ellipse)
[P_xy x0 y0] = projellisa(valP,[],'r');


end














