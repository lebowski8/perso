clear all
close all
home








% définition des variables
A1 = [0 1 ; -1 -2];
model1 = [A1 [0 ; 0] ; 0 0 1]

model2 = [-1 1 ; -1 -2];

model = {model1};

for numModel = 1:length(model)
    
    mod = model{numModel};
    A = mod(1:2,1:2);
    P = sdpvar(3,3, 'symmetric');

    % définition des contraintes
    F = P(1:2,1:2)>0;
    %F = F + ((-[eye(2) A']*[zeros(2) P ; P zeros(2)]*[eye(2) ; A]) > 0);
    F1 = A'*P(1:2,1:2) + P(1:2,1:2)*A;
    F = F + (F1<0);
    % définition des options du solver
    options = sdpsettings('verbose',1,'solver','sedumi');

    % résolution du problème par le solver
    sol = solvesdp(F,[],options);
    checkset(F);

    % on vérifie que P est définie positive
    valP = value(P(1:2,1:2));
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
              ind_error=[cpt1 cpt2] 
           end
        end
    end

    for cpt1=1:numel(x1)
        for cpt2=1:numel(x2)
           mydV=[x1(cpt1) x2(cpt2)]*A'*valP*[x1(cpt1);x2(cpt2)] + [x1(cpt1) x2(cpt2)]*valP*A*[x1(cpt1);x2(cpt2)];
           mydLyap(cpt1,cpt2) = mydV; 
        end
    end
    
%     x3=-10:1:10;
% 
%     for cpt1=1:numel(x1)
%         for cpt2=1:numel(x2)
%             for cpt3 = 1:numel(x3)
%                 myLyap(cpt1,cpt2,cpt3) = [x1(cpt1) x2(cpt2) x3(cpt3)]*valP*[x1(cpt1) ; x2(cpt2) ; x3(cpt3)];
%                 if myLyap(cpt1,cpt2,cpt3)<0
%                     ind_error=[cpt1 cpt2 cpt3] 
%                 end
%             end
%         end
%     end
% 
%     for cpt1=1:numel(x1)
%         for cpt2=1:numel(x2)
%             for cpt3 = 1:numel(x3)
%                mydV=[x1(cpt1) x2(cpt2) x3(cpt3)]*A'*valP*[x1(cpt1);x2(cpt2) ; x3(cpt3)] + [x1(cpt1) x2(cpt2) x3(cpt3)]*valP*A*[x1(cpt1);x2(cpt2) ; x3(cpt3)];
%                mydLyap(cpt1,cpt2,cpt3) = mydV; 
%             end
%         end
%     end


    %% Vérification
    if(all(all(myLyap>=0)))
        disp('V(x) est positive pour tout x')
    end
    
    if(all(all(mydLyap<=0)))
        disp('dV(x)/dx est négative pour tout x')
    end
    %% Affichage 



    % figure
    % h = meshc(x1,x2,myLyap);
    % title('V(x)')
    % xlabel('X(1)')
    % ylabel('X(2)')
    % zlabel('V(X)')


    figure(1)
    [C h] = contour(x1,x2,myLyap);
    title('projection de V(x)')
    hold on
    grid on
    xlabel('X(1)')
    ylabel('X(2)')

    % hold on
    % 
    % plot(C(1,2:34),C(2,2:34),'r')
    % plot(C(1,36:84),C(2,36:84),'g')
    % plot(C(1,86:146),C(2,86:146),'black')

end


