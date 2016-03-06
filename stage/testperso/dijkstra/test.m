clear all
close all
home


% Définition de la mission

Hmax = 10;
Lmax = 14;


zone1 = [2 2 6 6 ; 5 8 8 5]';
zone2 = [4 4 6.5 6.5 ; 0 2 2 0]';
zone3 = [7 7 8.5 8.5 8 8; 3 4 4 2 2 3]';
zone4 = [9.5 9.5 11.5 11.5 ; 0 2 2 0]';
%zone5 = [10 11 11 10 ; 10 10 9 9]';

zone = {zone1, zone2, zone3, zone4};
% zone = {zone1, zone2, zone3}

figure
hold on

fill(zone1(:,1),zone1(:,2),'r')
fill(zone2(:,1),zone2(:,2),'b')
fill(zone3(:,1),zone3(:,2),'g')
fill(zone4(:,1),zone4(:,2),'y')
%fill(zone5(:,1),zone5(:,2),'black')

axis([-1 Lmax+1 -1 Hmax+1])
xlabel('distance (km)')
ylabel('hauteur (km)')
grid on

text(0,-0.5,[' Départ']);
plot(0, 0,'r.','MarkerSize',20) % On affiche les points
text(12,4,['  Arrivé']);
plot(12, 4,'r.','MarkerSize',20) % On affiche les points



%% Création de la grille de point

% grille totale
for i = 1:(Lmax+1)
    for j = 1:(Hmax+1)
        n = ((i-1)*(Hmax+1) + j);
        nodes( n , :) = [n i-1 j-1]; 
    end
end

for z = 1:length(zone)
    zoneZ = zone{z};
    for i = 1:length(zoneZ)
    indx = find((zoneZ(i,1)==nodes(:,2)));

    indy = find((zoneZ(i,2)==nodes(:,3)));

    ind = intersect(indx,indy);

    nodes(ind,:) = [];
    end
end

% plot(nodes(:,2), nodes(:,3),'k.')
% 
% for n = 1:(Lmax*Hmax)
%         text(nodes(n,2),nodes(n,3),[' ' num2str(nodes(n,1))]); % on numérote les noeuds
% end
%%
% Suppression des points dans les zones


% suppression des noeuds à la main -> il faudra l'automatiser
% pour ne pas avoir de soucis avec les indices, il faut dabord faire
% tourner sans ce qui suit, puis utiliser les noeuds afficher pour
% supprimer les lignes à la main....
%del = {40 41 51 52 56 57 67 68 92 93 102 103 104 111 122 123 124}
% del = [37, 38, 45, 46, 49 50 58 59 60 80 81 90 91 92 99 109 110 111];
del = [29 30 39 40 41 42 46 50 51 52 53 56 57 58 61 62 63 64 67 68 69 73 74 93 111 112 113 122 123 124];
% for i = 1:length(del)
%    nodes( (del{i}-i+1),:) = []; 
% end

for i = 1:length(del)
    ind = find(nodes(:,1) == del(i));
    nodes(ind) = 0;
end

for i = 1:length(del)
    ind = find(nodes(:,1) == 0);
    nodes(ind,:) = [];
end


% plot(nodes(:,2), nodes(:,3),'k.') % On affiche les points

for n = 1:length(nodes)
%         text(nodes(n,2),nodes(n,3),[' ' num2str(nodes(n,1))]); % a
%         utiliser pour sélectionner les noeuds à supprimer
%         text(nodes(n,2),nodes(n,3),[' ' num2str(n)]); % a utiliser pour
%         selectionner les pointsde départs et arrivés après avoir
%         supprimés correctement les noeuds
end

nodes(:,1) = 1:length(nodes); % on met à jour la première colonne pour que les indices corresondent

%%
%nodes = [(1:10) ; 10*rand(2,10)]';

%segments = [(1:17) ; floor(1:0.5:9) ; ceil(2:0.5:10)]'

% segments = [1 1 13
%             2 13 25
%             3 25 35
%             4 35 44
%             5 44 53
%             6 44 54
%             7 53 61
%             8 61 70
%             9 70 78
%             10 78 89
%             11 89 100
%             12 100 112
%             13 112 122
%             14 122 121
%             15 121 131
%             16 54 63
%             17 63 71
%             18 71 82
%             19 82 93
%             20 93 102
%             21 102 112];
        
    indx = find( (nodes(62,2)+1) == nodes(:,2) );
    indy = find( (nodes(62,3)) == nodes(:,3) );
    ind = intersect(indx,indy);

k = 1;
for i = 1:length(nodes)
    ni = nodes(i,2:3); 
    xi = ni(1);
    yi = ni(2);

    if yi == 0 % Si on est à hauteur mini
        % recherche du noeud à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
        
        % recherche du noeud en haut à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)+1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end

    elseif yi == Hmax % Si on est à hauteur max
        % recherche du noeud à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
        
        % recherche du noeud en bas à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)-1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
        
        % recherche du noeud en bas
        indx = find( (nodes(i,2)) == nodes(:,2) );
        indy = find( (nodes(i,3)-1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
    else
        % recherche du noeud en haut à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)+1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end

        % recherche du noeud à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
        
        % recherche du noeud en bas à droite
        indx = find( (nodes(i,2)+1) == nodes(:,2) );
        indy = find( (nodes(i,3)-1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
        
        % recherche du noeud en bas 
        indx = find( (nodes(i,2)) == nodes(:,2) );
        indy = find( (nodes(i,3)-1) == nodes(:,3) );
        ind = intersect(indx,indy);
        if ~isempty(ind)
            segments(k,:) = [k i ind];
            k = k+1; 
        end
    end

end
        
        
        
%% Dijkstra
%figure 
%plot(nodes(:,2), nodes(:,3),'k.')
hold on

% for s = 1:length(segments)
% %     if (s <= 10) 
% %         text(nodes(s,2),nodes(s,3),[' ' num2str(s)]); 
% %     end
%     plot(nodes(segments(s,2:3)',2),nodes(segments(s,2:3)',3),'b');
% end
disp('Résultat Dijkstra')
 [d, p] = dijkstra(nodes, segments, 1, 97) % ici 1 est le noeud de départ et 124 celui d'arrivé

for n = 2:length(p)
    plot(nodes(p(n-1:n),2),nodes(p(n-1:n),3),'r-.','linewidth',2);
end
hold off;

%% Bellman-ford

%W = [.41 .99 .51 .32 .15 .45 .38 .32 .36 .29 .21];
%DG = sparse([6 1 2 2 3 4 4 5 5 6 1],[2 6 3 5 4 1 6 3 4 3 5],W);

DG = sparse(segments(:,2),segments(:,3),ones(1,length(segments)));

[m,n,p,D,tail,head,W] = Initialize(DG);
[p, D,iter] = BFMSpathOT(DG,1,97);

disp('Résultat Bellman-Ford')
if iter == 0
    disp('Pas de chemin')
elseif length(p) == length(D)
    chemin = [(1:length(p))' p D]
else
    p
    D
end

hold on
for n = 2:length(p)
    plot(nodes(p(n-1:n),2),nodes(p(n-1:n),3),'b-.','linewidth',2);
end



