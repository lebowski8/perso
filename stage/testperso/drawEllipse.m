function drawEllipse()
% Cette fonction prend en paramètre une matrice 2x2 et calcul les vecteurs
% propres / valeurs propres afin de tracer l'ellipse associée

A = [0.6038 0.5637 ; 0.5637 1.7960]

[V D] = eig(A);

S = sortrows([V' diag(D)],3);

Vmin = S(1,1:2)';
lmin = S(1,3);

Vmax = S(2,1:2)';
lmax = S(2,3);

%[X Y] = calculateEllipse(0, 0, lmin/2,lmax/2,-atan2(Vmin(2),Vmin(1)));



% figure
% grid on
% axis equal

% plot(X,Y)
% hold on
% plot([0 Vmin(1)],[0 Vmin(2)],'g')
% plot([0 Vmax(1)],[0 Vmax(2)])
% 
% x = linspace(-2,2,100);
% plot(x,sqrt((1-lmin.*x.^2)/lmax))

end