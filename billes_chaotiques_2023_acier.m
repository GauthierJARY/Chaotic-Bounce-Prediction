function [ Gamma, vphi ] = billes_chaotiques_2023_bois(name, f)


%% Données initiales
g = 9.81;
N=200;   % Paramètre de lissage pour l'interpolation, nombre de points par période
T=1/f;
omega=2*pi/T;
npp=8; % Nombre de phase par groupe
Np = N*npp; % nombre de points par groupe de phase

%% Chargement des données d'acquisition
a = load(name);
t=a(:,1);
zb=a(:,3);
zp=a(:,2);
facq=1/mean(diff(t));

%% Recentrage des données et affichage
zb=(-zb+mean(zb))+1.0474; % Valeur ajoutée à la main pour une bonne lecture graphique     % FILL HERE
zp=zp-mean(zp);

figure(1); 
plot(t,zb,'r-'); hold on;
plot(t,zp,'b-');
hold off ;


%% Interpolation des données et Moyenne de phase

ti=0:T/N:max(t);                                              % La nouvelle fréquence de hachage du signal est de N/T
zbi=interp1(t,zb,ti);
zpi=interp1(t,zp,ti);
%plot(ti,zbi,'g-'); 
%plot(ti,zpi,'y-');

% Nombre de période pour l'affichage de la moyenne de phase
np=floor(max(t)/(T*npp));

% Calcul de la moyenne de phase
zbp=zbi(1:Np);
zpp=zpi(1:Np);
for i=2:np
    zbp=zbp+zbi((i-1)*Np+1:i*Np); 
    zpp=zpp+zpi((i-1)*Np+1:i*Np);
end

% Normalisation de la moyenne
zbp=zbp/np;
zpp=zpp/np;

% Mise en forme graphique
zbp = zbp - (min(zbp) - max(zpp));

% Affichage de la moyenne de phase
tp=ti(1:Np);
figure(2);
plot(tp,zbp,'r-'); hold on;
plot(tp,zpp,'g-');

% Calibrage de l'amplitude, du temps et conversion en m
A = std(zpp)*sqrt(2); % On calcule la variance pour utiliser tous les points
[zppmax,id] = max(zpp);
to = tp(id);

    % ********************************************************************
    % De par nos observations, on remarque qu'une amplitude de 0.2V
    % correspond à 0.2mm, On a donc A en m sans "convertir" directement
    % ********************************************************************

% Fit du plateau
f = ezfit(tp,zpp,['z(t) = c + A*sin(',num2str(omega), '*(t-to)); c=0;A=',num2str(A),'; to=',num2str(to)]);
plot(tp,f.m(2)+f.m(1)*sin(omega*(tp-f.m(3))),'k-');


%% Extraction des rebonds 

zbps=smooth(zbp,60,'lowess');
plot(tp,zbps,'ko','MarkerSize',3); hold on;
a = imregionalmin(zbps)  % Crée une matrice de la taille de zbp avec 1 quand il y a un min, 0 sinon
c = imregionalmax(zbps);

%FILTRE
% for i=1:length(a)   
%     if a(i) > mean(zbp) && a(i)==1
%         a(i)=0;
%     end
% end
% a

ind = find(a==1);
ind = ind(2:end-1);
cm = mean(zbps(c));

% On affiche les points minimums
plot(tp(ind),zbp(ind),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');

% Définir les intervalles autour des minima
    % Ce choix est subjectif, car la dérivée seconde d'un signal
    % expérmiental serait extrêmement bruité

%% test pour un point / boucle for
ne=20; % Nombre de points exclus
G = -g/2*1e2; % Accélération de la pesanteur en terme d'amplitude de mouvement

% On regarde le premier minimum

% Définition de l'intervalle
interg = [max(ind(1)-N+ne,0+ne):ind(1)-ne];
interd = [ind(1)+ne:ind(2)-ne];

tmg=tp(ind(1))-T/2;
tmd=tp(ind(1))+T/2;

plot(tp(interg),zbp(interg),'ko','MarkerSize',8);
plot(tp(interd),zbp(interd),'ko','MarkerSize',8);

% Fit parabolique sur chacun des intervalles
fg = ezfit(tp(interg),zbp(interg),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmg),';c=',num2str(cm)]);
fd = ezfit(tp(interd),zbp(interd),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmd),';c=',num2str(cm)]);

% Trouver le point d'intersection
alpha = -2*G*(fd.m(2)-fg.m(2));
beta = fg.m(1) - fd.m(1) + G*(fg.m(2)^2-fd.m(2)^2);

tc = beta/alpha;

% Affichage de l'intersection
plot(tp,G*(tp-fg.m(2)).^2 + fg.m(1),'b-'); hold on;
plot(tp,G*(tp-fd.m(2)).^2 + fd.m(1),'m-');
plot(tc(1),G*(tc(1)-fd.m(2)).^2 + fd.m(1),'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k');

for i=2:length(ind)-1
    
    % Définition de l'intervalle
    interg = [ind(i-1)+ne:ind(i)-ne];
    interd = [ind(i)+ne:ind(i+1)-ne];
    
    tmg=tp(ind(i))-T/2;
    tmd=tp(ind(i))+T/2;
    
   plot(tp(interg),zbp(interg),'ko','MarkerSize',8);
   plot(tp(interd),zbp(interd),'ko','MarkerSize',8);
    
    % Fit parabolique sur chacun des intervalles
    fg = ezfit(tp(interg),zbp(interg),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmg),';c=',num2str(cm)]);
    fd = ezfit(tp(interd),zbp(interd),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmd),';c=',num2str(cm)]);
    
    % Trouver le point d'intersection
    alpha = -2*G*(fd.m(2)-fg.m(2));
    beta = fg.m(1) - fd.m(1) + G*(fg.m(2)^2-fd.m(2)^2);
    
    tc(i) = beta/alpha;
    
    plot(tp,G*(tp-fg.m(2)).^2 + fg.m(1),'b-'); hold on;
    plot(tp,G*(tp-fd.m(2)).^2 + fd.m(1),'m-');
    plot(tc(i),G*(tc(i)-fd.m(2)).^2 + fd.m(1),'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k');
    
end
  
%On regarde le dernier minimum

% Définition de l'intervalle
imax = length(ind) ;
interg = [ind(imax-1)+ne:ind(imax)-ne];
interd = [ind(imax)+ne:min(ind(imax)+N-ne,length(tp)-ne)];

tmg=tp(ind(imax))-T/2;
tmd=tp(ind(imax))+T/2;

% Fit parabolique sur chacun des intervalles
fg = ezfit(tp(interg),zbp(interg),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmg),';c=',num2str(cm)]);
fd = ezfit(tp(interd),zbp(interd),['z(t) = ',num2str(G),'*(t-t0)^2 + c; t0=',num2str(tmd),';c=',num2str(cm)]);

% Trouver le point d'intersection
alpha = -2*G*(fd.m(2)-fg.m(2));
beta = fg.m(1) - fd.m(1) + G*(fg.m(2)^2-fd.m(2)^2);

tc(imax) = beta/alpha;

plot(tp,G*(tp-fg.m(2)).^2 + fg.m(1),'b-'); hold on;
plot(tp,G*(tp-fd.m(2)).^2 + fd.m(1),'m-');
plot(tc(imax),G*(tc(imax)-fd.m(2)).^2 + fd.m(1),'o','MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k');

figure(2); ylim([-1 2.5]);

t0 = f.m(3);
vphi = omega*(tc-t0);
%vphi = mod(vphi,2*pi);
vphi = wrapTo2Pi(vphi) ;

% Accélération adimensionnée, paramètre de notre problème
Gamma = (f.m(1)*0.01*omega^2)/9.81 ;
hold off ;
%     % Affichage des fits
%     plot(tp(interg),G*(tp(interg)-fg.m(2)).^2 + fg.m(1),'b-');
%     plot(tp(interd),G*(tp(interd)-fd.m(2)).^2 + fd.m(1),'m-');
        
%     % Affichage du point d'intersection
%     figure(3)
%     plot(tp,G*(tp-fg.m(2)).^2 + fg.m(1),'b-'); hold on;
%     plot(tp,G*(tp-fd.m(2)).^2 + fd.m(1),'m-');
%     plot(tc,G*(tc-fd.m(2)).^2 + fd.m(1),'o','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k');


end

