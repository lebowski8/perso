function [valP] = calcMatriceLyap(hold,descent,climb)

Ahold = hold.a;
Bhold = hold.b;

Adescent = descent.a;
Bdescent = descent.b;

Aclimb = climb.a;
Bclimb = climb.b;

% P00 : matrice pour switch de maintien vers maintien
% P01 : matrice pour switch de maintien vers descent
% P02 : matrice pour switch de maintien vers climb

P00 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P01 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P02 = sdpvar(length(Ahold),length(Ahold), 'symmetric');

P11 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P10 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P12 = sdpvar(length(Ahold),length(Ahold), 'symmetric');

P22 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P20 = sdpvar(length(Ahold),length(Ahold), 'symmetric');
P21 = sdpvar(length(Ahold),length(Ahold), 'symmetric');

% toutes les matrices doivent être définie positive
    F = P00>0;
    F = F + (P01>0);
    F = F + (P02>0);
    
    F = F + (P11>0);
    F = F + (P10>0);
    F = F + (P12>0);

    F = F + (P22>0);
    F = F + (P20>0);
    F = F + (P21>0);

% toutes les dérivées doivent être définies négative

    F = F + ([(Ahold'*P00*Ahold-P00) , (Ahold'*P00*Bhold) ; (Bhold'*P00*Ahold) , (Bhold'*P00*Bhold-eye(size(Bhold,2)))]<0);   
    F = F + ([(Ahold'*P01*Ahold-P01) , (Ahold'*P01*Bhold) ; (Bhold'*P01*Ahold) , (Bhold'*P01*Bhold-eye(size(Bhold,2)))]<0);   
    F = F + ([(Ahold'*P02*Ahold-P02) , (Ahold'*P02*Bhold) ; (Bhold'*P02*Ahold) , (Bhold'*P02*Bhold-eye(size(Bhold,2)))]<0);   
    
    F = F + ([(Adescent'*P11*Adescent-P11) , (Adescent'*P11*Bdescent) ; (Bdescent'*P11*Adescent) , (Bdescent'*P11*Bdescent-eye(size(Bdescent,2)))]<0);   
    F = F + ([(Adescent'*P10*Adescent-P10) , (Adescent'*P10*Bdescent) ; (Bdescent'*P10*Adescent) , (Bdescent'*P10*Bdescent-eye(size(Bdescent,2)))]<0);   
    F = F + ([(Adescent'*P12*Adescent-P12) , (Adescent'*P12*Bdescent) ; (Bdescent'*P12*Adescent) , (Bdescent'*P12*Bdescent-eye(size(Bdescent,2)))]<0);   

    F = F + ([(Aclimb'*P22*Aclimb-P22) , (Aclimb'*P22*Bclimb) ; (Bclimb'*P22*Aclimb) , (Bclimb'*P22*Bclimb-eye(size(Bclimb,2)))]<0);   
    F = F + ([(Aclimb'*P20*Aclimb-P20) , (Aclimb'*P20*Bclimb) ; (Bclimb'*P20*Aclimb) , (Bclimb'*P20*Bclimb-eye(size(Bclimb,2)))]<0);   
    F = F + ([(Aclimb'*P21*Aclimb-P21) , (Aclimb'*P21*Bclimb) ; (Bclimb'*P21*Aclimb) , (Bclimb'*P21*Bclimb-eye(size(Bclimb,2)))]<0);   
   
% entre chaque switch la fonction V doit décroitre, donc les matrices P
% doivent être plus "petite" : P00 - P01 semi-défine positive
    F = F + (P00 >= P01);
        F = F + (P01 >= P11);
            F = F + (P11 >= P10);
            F = F + (P11 >= P12);

        F = F + (P01 >= P10);
            F = F + (P10 >= P00);
            F = F + (P10 >= P01);
            F = F + (P10 >= P02);

        F = F + (P01 >= P12);
            F = F + (P12 >= P20);
            F = F + (P12 >= P21);
            F = F + (P12 >= P22);

    F = F + (P00 >= P02);
        F = F + (P02 >= P20);
            F = F + (P20 >= P00);
            F = F + (P20 >= P01);
            F = F + (P20 >= P02);
            
        F = F + (P02 >= P21);
            F = F + (P21 >= P10);
            F = F + (P21 >= P11);
            F = F + (P21 >= P12);
            
        F = F + (P02 >= P22);
            F = F + (P22 >= P20);
            F = F + (P22 >= P21);
                



% définition des options du solver
options = sdpsettings('verbose',1,'solver','lmilab');

% résolution du problème par le solver
sol = solvesdp(F,[],options);
checkset(F);

valP = {value(P00), value(P01), value(P02), value(P11), value(P10), value(P12), value(P22), value(P20), value(P21)};

%% on vérifie que P est définie positive (donc V positive pour tout x)
% valP = value(P);
% spectreP = eig(valP);
% flagP = 0;
% for i = 1:rank(valP)
% 	if spectreP(i) <= 0 
%         disp(spectreP(i))
%         flagP = 1;
%     end
% end
% 
% if flagP == 1
% 	disp('P : non définie positive')
% else
% 	disp('P : définie positive')
% end


end














