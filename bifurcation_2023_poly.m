clear all ; close all ; clc ;

% Nomenclature : test030-001.txt ,
% 030 pour la fréquence de 30Hz
% 001 pour le premier échantillon


list=rdir(['serie0703_30_*.txt']);

figure(10); hold off;
figure  
% periode triple pour 20 à 24
% il faut reprendre les derniers avec des bornes plus petites

% faire 2 scripts avec poly triple et poly double
for i=20:24 
    name = list{i}(20:39);
    strf = name(3:4) ;
    f = str2double(strf);
    f=30;
    [ Gamma, vphi ] = billes_chaotiques_2023_poly(name, f) ;
pause;
    figure(10); hold on; 
    plot(Gamma,vphi,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6);
end
%%
alpha=0.815
Gam=0:0.01:2;
phi_theo=acos((1-alpha)/(1+alpha)*pi./Gam)
plot(Gam,phi_theo,'k-');
%%
%alpha=0.94
Gamma1= pi *   (        ( (1-alpha)/(1+alpha) )^2    +    ( 2*(1+alpha^2)/(pi*( (1+alpha)^2 ) ) )^2        )^0.5
Gamma2=pi*(1+4/(pi^2)*((1+alpha^2)/(1-alpha^2))^2)^0.5*((1-alpha)/(1+alpha))
plot([Gamma2 Gamma2],[0 2.5],'k--')

%%
title('Courbe de bifurcation avec une balle de polystyrène à 30 Hz');
xlabel('Gamma');
ylabel('Amplitude');