% Analise de tempo computacional em funcao do numero de agentes

t = [4.7748;
    6.7888;
    8.7236;
    9.6433
    11.1707;
    14.4767;
    16.8881
   18.3683;
   23.4542;
   32.2390;
   55.5227];
ag = [3;4;5;6;7;8;9;10;12;15;20];

figure (15)
hold on
grid on;
% xlim([0 t(end)])
plot(ag,t,'-o','LineWidth',1.5); hold on;
hold on;  box on;
ylabel('Tempo [s]'); xlabel('Numero de agentes')