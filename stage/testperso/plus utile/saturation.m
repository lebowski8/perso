clear;

%systeme
Ap=0.1;
Bp=1;
Cp=1;

l=length(Ap);
%PI
%Ar=0;
%Br=0.2;
%Cr=1;
%Dr=2;

%FORE
Ar=-1;
Br=1;
Cr=1;
Dr=0;
m=length(Ar);
nc=length(Ar);
%syst complet

A=[Ap-Bp*Dr*Cp Bp*Cr;
    -Br*Cp Ar];
n=length(A);
B=[Bp;
    0];
K=[-Dr*Cp Cr];
Aj=[eye(l) zeros(l,m);
    zeros(m,l) zeros(m)];
u0=1;


M1=[0 1; 1 0];
Q=[-Cp zeros(nc,1);
    -Dr*Cp Cr];

M=Q'*M1*Q;



P=sdpvar(n,n,'sym');
G=sdpvar(1,n);
% Y2=sdpvar(1,n);
S=sdpvar(1,1);

tauf=sdpvar(1);
taur=sdpvar(1);



lmic=[A'*P+P*A+tauf*M P*B-G';
    B'*P-G -2*S];

lmid=Aj'*P*Aj-P-taur*M;

lmi1=[P K'*S-G';
    S*K-G u0^2];

%  lmi1=[W W*K'-(Y1)';
%      K*W-(Y1) u0^2];
%condition initiale
delta=sdpvar(1,1);
v{1}=[zeros(m-1,1);1;0];
v{2}=[zeros(m-1,1);-1;10];
v{3}=[zeros(m-1,1);-1;0];
v{4}=[zeros(m-1,1);1;-10];

for i=[1 3]
    opt{i}=[delta-v{i}'*P*v{i}];
end
%contraintes
cont=(P>0)+(S>0)+(lmic<0)+(lmi1>=0)+(lmid<1e-9*eye(length(lmid)));
cont=cont+(tauf>0)+(taur>0);
cont=cont+(delta>0);
for i=[1]
  cont=cont+(opt{i}>=0);
end

%resolution
option = sdpsettings('solver','sdpt3');
solvesdp(cont,delta-4*S,option);

[p d]=checkset(cont);
if min(p)<0  
  disp('erreur');
end

t=double(P);
s=inv(double(S));
figure(30);hold on
projellisa(t*s^2,[],'g');


%deuxièmeresolution
% r=[1 0];
%  P1=sdpvar(n,n,'sym');
%  a1=sdpvar(1);
%  lmi2=Aj'*(t*s^2)*Aj-P1;
%  lmi3=Aj'*(t*s^2)*Aj-a1*r'*r
%  pb=set(a1>0)+set(lmi3(1,1)<=0);
%  option = sdpsettings('solver','lmilab');
%  solvesdp(pb,a1,option);
%  
%  [p d]=checkset(cont);
%  if min(p)<0  
%    disp('erreur');
%  end
% % 
%   eli=1/sqrt(double(a1))
%  % figure(10);hold on
%   %projellisa(eli,'--','b');
