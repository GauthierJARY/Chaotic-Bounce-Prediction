a=load('serie_020.txt');

t=a(:,1);
zp=a(:,2);
zb=-a(:,3);

figure(1); hold off; 
plot(t,zp,'b.');
hold on;
plot(t,zb,'r-');


%%