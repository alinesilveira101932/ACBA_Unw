%% Comparação horizonte de prefição
% Com e sem tratamento de unwinding
% 05/03/2026
% Aline Isabel
% ----------------------------------------------
clc;clear;close all;
Restricoes = [35; 30; 25; 20; 15; 10; 5; 3; 2];
minHComTratUnwind = [2; 2;2; 2; 2; 2; 4; 17 ; 34];
minHSemTratUnwind = [2; 2;2; 4; 6; 12; 30; 61; 105];

figure (1)
hold on
grid on; hold on;

plot(Restricoes,minHComTratUnwind,'-or','LineWidth',2)
plot(Restricoes,minHSemTratUnwind,'-ob','LineWidth',2)
xlabel('u_{max}')
ylabel('H_{min}')
legend('ComTratamento', 'SemTratamento')
hold on;
% plot([0 kend],[umax1 umax1],'--k','linewidth',1,'handlevisibility','off')
% text(0.1,0.18,'Limitante de x2','Color','r','FontWeight','bold')
