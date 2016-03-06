% function projellisa(P,linetype,couleur,holdflag)
%
% visualisation de l'ellipsoide {x'*P*x <= 1} dans toutes les projections

function [P_xy x_cor y_cor]=projellisa(P,linetype,couleur,holdflag)

if nargin == 1,
 linetype = '-';
 couleur = 'k';
end;
if nargin < 4,
 holdflag = 0;
end;

%keyboard
n = size(P,1);

if n == 2,
 nrow = 1; ncol = 1;
elseif n == 3,
 nrow = 3; ncol = 1;
elseif n == 4,
 nrow = 3; ncol = 2;
elseif n == 5,
 nrow = 5; ncol = 2;
elseif n == 6,
 nrow = 5; ncol = 3;
elseif n == 7,
 nrow = 7; ncol = 3;
elseif n == 10
    nrow = 9 ; ncol = 5; 
elseif n == 11
    nrow = 11 ; ncol = 5; 
else
 error('invalid dimension for P');
end;

x1 = 1; x2 = 2;
 %figure
for i = 1:n*(n-1)/2
%  hold on;
%  grid on
%  axis equal
%  subplot(nrow, ncol, i);




 % -------------------------
 % Modifs pour la projection
 ind1 = [x1 x2];
 ind2 = [1:x1-1 x1+1:x2-1 x2+1:n];
 A = P(ind2, ind2);
 B = P(ind1, ind2);
 C = P(ind1, ind1);
%  [x_cor y_cor]=plotellisa(C-B*inv(A)*B',zeros(n,1),linetype,couleur);
 
 % affichage romain
 if ( (x1==4) && (x2==8) )   % on ne cherche que les ellipses qui nous interessent
     P_48 = C-B*inv(A)*B';
      [x_48 y_48]=plotellisa(P_48,zeros(n,1),linetype,couleur);

%      hold on
%      grid on
%      %axis equal
%      subplot(3, 1, 1);
%      plot(x_cor,y_cor)
%      xlabel('Va')
%      ylabel('h')
     
 elseif ( (x1==4) && (x2==5) )
     P_45 = C-B*inv(A)*B';
       [x_45 y_45]=plotellisa(P_45,zeros(n,1),linetype,couleur);
% 
%     hold on
%     grid on
%     %axis equal
%     subplot(3, 1, 2);
%     plot(x_cor,y_cor)
%     xlabel('Va')
%     ylabel('alpha')
    
 elseif ( (x1==5) && (x2==8) )
     P_58 = C-B*inv(A)*B';
      [x_58 y_58]=plotellisa(P_58,zeros(n,1),linetype,couleur);

%     hold on
%     grid on
%     %axis equal
%     subplot(3, 1, 3);
%     plot(x_cor*0.005,y_cor)
%     xlabel('alpha')
%     ylabel('h')
    
 end

 %disp(['x' int2str(x1) '/' ['x' int2str(x2)]]);
 
% grid on
% axis equal
%  xlabel(['x' int2str(x1)]);
%  ylabel(['x' int2str(x2)]);
 % -------------------------
%  xlabel(['x' int2str(x1)]);
%  ylabel(['x' int2str(x2)]);
%  grid on
%  axis equal
 if x2 < n,
  x2 = x2 + 1;
 else
  x1 = x1 + 1;
  x2 = x1 + 1;
 end;  

end;
 P_xy = {P_48 ; P_45 ; P_58};
 x_cor = [x_48 ; x_45 ; x_58];
 y_cor = [y_48 ; y_45 ; y_58];
end

