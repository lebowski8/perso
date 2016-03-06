% Test chauffage
close all
clear all
Ts=0.005;                % Pas de temps
T0=18;                  % Temperature initiale
Tlim=22;              % Température à ne pas dépasser en maxi
Tmini = 15;             % Température à ne pas dépasser en mini
tsim=5;                % Temps de simulation
%
t = 0:Ts:tsim-Ts;
in=(Tlim+1)*ones(tsim/Ts,1);  % Entree positive
X=zeros(tsim/Ts,1);           % Vecteur d'état
Y=zeros(tsim/Ts,1);           % Vecteur de sortie
X(1)=T0;                      % Etat initial (pris comme la sortie initiale)

%%

for k=1:tsim/Ts
   if T0<=Tlim                % Système avec entrée (température non atteinte)
       A=-1;B=1;C=1;D=0;       % Représentation d'état (continue)
       sys1=ss(A,B,C,D);
       sys1=c2d(sys1,Ts);      % Passage en discret
       X(k+1)=sys1.A*X(k)+sys1.B*in(k); 
       Y(k)=sys1.C*X(k)+sys1.D*in(k);
       T0=Y(k);
   else
       A=-1;B=0;C=1;D=0;        % Système sans entrée (switch si T>Tlim)
       sys2=ss(A,B,C,D);
       sys2=c2d(sys2,Ts);
       X(k+1)=sys2.A*X(k)+sys2.B*in(k);
       Y(k)=sys2.C*X(k)+sys2.D*in(k);
       if Y(k) <= Tmini
            T0=Y(k);
       end
   end
end

%%

figure
plot(t',Y)
grid on
xlabel('temps')
ylabel('temperature')