clc; clear; close all;

%% Parametros

Parametros_yalmip

% Objeto
robo = robot_model;
robo.A = A;
robo.B = B;
robo.x0 = [pi/3;10;pi/10;0;pi/3;-5];
robo.State_ref = Nx;
robo.Input_ref = Nu;
robo.Term_bo = bo;
robo.Term_Or = Or;
robo.Term_Sf = Sf;

%% Otimizacao
maxN = N;

cenario = probRestricoes;
cenario.Robot = robo;
cenario.rho = R_barra;
cenario.rho_x = Q_barra;
cenario.P = P_barra;
cenario.maxN = maxN;
cenario.create_yalmip_prob;


%% Otimizando
kend = 100;

for i = 1:kend

cenario.solve;
cenario.SolTime;
xk(:,i) = cenario.Robot.X(:,1);
uk(:,i) = cenario.Robot.U(:,1);
rbark(:,i) = cenario.Robot.rbar_var;
Nbark(:,i) = cenario.Robot.N_var;
Rbar(:,i) = [rbark(:,i); rbark(:,i) + 2*pi*Nbark(1,i); rbark(:,i) + 2*pi*Nbark(2,i)];
robo.evolve_state(uk(:,i));
end


%% Plots
x = [xk cenario.Robot.X(:,1)];
t = 0.1*[0:kend];

figure (104)
subplot(4,1,1)
grid on; 
plot([0 t(end)],[0 0],'LineStyle',':','Color',[0.39,0.83,0.07],'LineWidth',2)
xlim([0 t(end)])
xlabel('Tempo, s');  box on;
legend('N_1')

subplot(4,1,2)
stairs(t(2:end),Nbark(1,:),'LineStyle',':','Color',[0.07,0.62,1.00],'LineWidth',2); hold on;
xlim([0 t(end)]); grid on; box on;
xlabel('Tempo, s');  
legend('N_2')

subplot(4,1,3)
stairs(t(2:end),Nbark(2,:),'LineStyle',':','Color',[0.93,0.69,0.13],'LineWidth',2); hold on;grid on;
xlim([0 t(end)])
xlabel('Tempo, s');   box on;
legend('N_3')

subplot(4,1,4)
plot(t(2:end),rbark,':k','LineWidth',2); hold on;grid on;
xlim([0 t(end)])
xlabel('Tempo, s');   box on;
legend('$\bar{r}$', 'Interpreter', 'latex')

check1 = abs(Rbar(1,end)- Rbar(2,end))/(2*pi)
check2 = abs(Rbar(1,end)- Rbar(3,end))/(2*pi)
check3 = abs(Rbar(3,end)- Rbar(2,end))/(2*pi)

figure (2)
hold on
grid on; xlim([0 t(end)])
plot(t(2:end),uk(1,:),'LineStyle',':','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t(2:end),uk(2,:),'LineStyle',':','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t(2:end),uk(3,:),'LineStyle',':','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u')

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x(1,:),'LineStyle',':','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(2,:),'LineStyle',':','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
hold on; box on;
stairs(t(2:end),Rbar(1,:),':k','linewidth',1)
xlabel('Tempo, s')
ylabel('Agente 1')
legend('$x_1$','$x_2$','$\bar{R}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')

figure (7)
hold on
grid on; 
plot(t,x(3,:),'LineStyle',':','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(4,:),'LineStyle',':','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),Rbar(2,:),':k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 2')
% legend('x_1','x_2')

figure (8)
hold on
grid on; 
plot(t,x(5,:),'LineStyle',':','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(6,:),'LineStyle',':','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),Rbar(3,:),':k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 3')


e12 = abs(exp(1j*x(1,:))-exp(1j*x(3,:))); e13 = abs(exp(1j*x(1,:))-exp(1j*x(5,:))); e23 = abs(exp(1j*x(3,:))-exp(1j*x(5,:)));

figure (13)
subplot(3,1,1)
hold on
grid on; xlim([0 t(end)])
plot(t,e12,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{12}')

subplot(3,1,2)
hold on
grid on; xlim([0 t(end)])
plot(t,e13,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{13}')
subplot(3,1,3)
hold on
grid on; xlim([0 t(end)])
plot(t,e23,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{23}')