
clc
format shortG

X = [0 0 0 0 0 0 0 0 0 0 0]'
u = [235 20]';
for i = 1:2
    Xnext=sys_hold.A*X+sys_hold.B*u;
    Y=sys_hold.C*X+sys_hold.D*u;

    X = Xnext
end

Xnext;
Y

%%
clc
format shortG

X = [0 0 0 235 0 0 0 1000 0 0]';
u = [235 5000]';
for i = 1:1
    Xnext=sys_climb.A*X+sys_climb.B*u;
    Y=sys_climb.C*X+sys_climb.D*u;

    X = Xnext;
end

Xnext
Y

%% 
%clc

X = [0 0 0 235 0 0 0 1500 0 0]';
% X = [18.6 0 0 64.7 -4.1 0 0 1077.3 60.1 0]';
u = [235 1500]';
for i = 1:1
    Xnext=sys_climbk5.A*X+sys_climbk5.B*u;
    Y=sys_climbk5.C*X+sys_climbk5.D*u;

    X = Xnext;
end

Xnext
Y

%% 
clc
clear X,Y;

X = [0 0 0 0 0 0 0 0 0 0]'
u = [235 20]';
for i = 1:5
    Y=sys_descentk5.C*X+sys_descentk5.D*u
    Xnext=sys_descentk5.A*X+sys_descentk5.B*u
    

    X = Xnext;
end

Xnext;
Y;

%% k10
clc
clear X,Y;
X = [0 0 0 0 0 0 0 0 0 0 0]'
u = [235 20]';



for i = 1:2
    Y=sys_holdk5.C*X+sys_holdk5.D*u
    Xnext=sys_holdk5.A*X+sys_holdk5.B*u    

    X = Xnext;
end

Xnext;
Y;

%%
clc
format shortG
u = [235 1499.9]';
% X = X(1:10);
X = [0 0 0 235 0 0 0 1000 0 0]'
for i = 1:4
    Xnext=sys_climbk5.A*X+sys_climbk5.B*u;
    Y=sys_climbk5.C*X+sys_climbk5.D*u

    X = Xnext;
end

u = [235 1125]';
for i = 1:1
    Xnext=sys_climbk5.A*X+sys_climbk5.B*u;
    Y=sys_climbk5.C*X+sys_climbk5.D*u

    X = Xnext;
end

Xnext;
Y;