%% MPC - Usando YALMIP - Tecnica Artigo Li2018 
% Autor: Aline Isabel
% Data: 05/2026
clc;clear;
% close all;
%--------------------------------------------------------------------
% Load parametros
Parametros_yalmip
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45

kend = 150;
 
% Estados iniciais
x1_0 = [pi/3;10];
x2_0 = [pi/10;0];
x3_0 = [pi/3;-5];
% 
% x1_0 = [0;0];
% x2_0 = [pi/5;0];
% x3_0 = [pi/2;0];

x(:,1) = [x1_0;x2_0;x3_0];

controller = mpc_controller(A,B,N,Q,R,S,na,nx,nu,Au,bu,I_na,beta); % Cria o problema

for k = 1:kend
    input = x(:,k);

    % Solucao
    [sol, diagnostic] = controller(input);
    u_sol = sol(N+2:end);
    x_sol = sol(1:N+1);
   
    % Opt
     u_opt(:,k) = u_sol{1};

    % Simulacao e evolucao dos estados
    
    % xini = x(:,k);
    % [t, xd] = ode45(@(t,x) simulation(x, u_opt(:,k),kron(I_na,Ac),kron(I_na,Bc)),[0 T], xini, options2); 
    % x(:,k+1) = xd(length(t),:)';
    x(:,k+1) = kron(I_na,A)*x(:,k)+kron(I_na,B)*u_opt(:,k);
   
end
% 
% function xdot = simulation(x,u,Ac,Bc)
% xdot = Ac*x+Bc*u;
% end

% tempo
t = 0.1*[0:kend];

%% Plots
figure (102)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),u_opt(1,:),'LineStyle',':','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t(1:end-1),u_opt(2,:),'LineStyle',':','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t(1:end-1),u_opt(3,:),'LineStyle',':','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u')
legend('Ag1', 'Ag2', 'Ag3')

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x(1,:),'LineStyle',':','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(2,:),'LineStyle',':','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
hold on; box on;
xlabel('Tempo, s')
ylabel('Agente 1')


figure (7)
hold on
grid on; 
plot(t,x(3,:),'LineStyle',':','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(4,:),'LineStyle',':','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
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
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 3')

figure (100)
hold on
grid on; xlim([0 t(end)])
plot(t,x(1,:),'LineStyle',':','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t,x(3,:),'LineStyle',':','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t,x(5,:),'LineStyle',':','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on; legend('Ag1', 'Ag2', 'Ag3')
xlabel('Tempo, s'); ylabel('x1')

figure (101)
hold on
grid on; xlim([0 t(end)])
plot(t,x(2,:),'LineStyle',':','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t,x(4,:),'LineStyle',':','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t,x(6,:),'LineStyle',':','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on; legend('Ag1', 'Ag2', 'Ag3')
xlabel('Tempo, s'); ylabel('x2')

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

f12 = abs(exp(1j*x(2,:))-exp(1j*x(4,:))); f13 = abs(exp(1j*x(2,:))-exp(1j*x(6,:))); f23 = abs(exp(1j*x(4,:))-exp(1j*x(6,:)));
figure (14)
subplot(3,1,1)
hold on
grid on; xlim([0 t(end)])
plot(t,f12,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('f_{12}')

subplot(3,1,2)
hold on
grid on; xlim([0 t(end)])
plot(t,f13,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('f_{13}')
subplot(3,1,3)
hold on
grid on; xlim([0 t(end)])
plot(t,f23,':b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('f_{23}')