clear all ; close all ; clc ;

% Nomenclature : test030-001.txt ,
% 030 pour la fr�quence de 30Hz
% 001 pour le premier �chantillon


list=rdir(['serie2102_20_*.txt']);

figure(10); hold off;
figure  
for i=1:length(list)
    name = list{i}(20:39);
    strf = name(3:4) ;
    f = str2double(strf);
    f=30;
    [ Gamma, vphi ] = billes_chaotiques_2023_acier(name, f) ;
pause;
    figure(10); hold on; 
    plot(Gamma,vphi,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6);
end
%%
alpha=0.944
Gam=0:0.01:2;
phi_theo=acos((1-alpha)/(1+alpha)*pi./Gam)
plot(Gam,phi_theo,'k-');
%%
alpha=0.94
Gamma1= pi *   (        ( (1-alpha)/(1+alpha) )^2    +    ( 2*(1+alpha^2)/(pi*( (1+alpha)^2 ) ) )^2        )^0.5
Gamma2=pi*(1+4/(pi^2)*((1+alpha^2)/(1-alpha^2))^2)^0.5*((1-alpha)/(1+alpha))
plot([Gamma2 Gamma2],[0 2.5],'k--')

%%
title('Courbe de rebonds avec une balle de Ping Pong � 30 Hz');
xlabel('Gamma');
ylabel('Amplitude');