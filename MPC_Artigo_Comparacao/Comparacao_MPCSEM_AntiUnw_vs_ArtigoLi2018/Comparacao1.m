%% MPC - Comparacao sem anti Unwinding  vs artigo li2018
% % Modelo artigo Li2018 - 5 estados e 2 entradas
% Autor: Aline Isabel
% Data: 24/07/2026
clc;clear;
% close all;
%--------------------------------------------------------------------
% Load parametros
Parametros_comparacao1
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45

kend = 150;
 
% Estados iniciais
for i = 1:na
x{i}(:,1) = xini(5*i-4:5*i);
end
rbar_old = [0;0];

controller = mpc_controller_comparacao(A,B,N,Q,R,P,W,Nx,Nu,bo,Or,Sf,na,nx,nu,Su,bu); % Cria o problema

for k = 1:kend
    input = [x{1}(:,k); x{2}(:,k); x{3}(:,k); x{4}(:,k); x{5}(:,k); rbar_old];

    % Solucao
    [sol, diagnostics] = controller(input);

    rbar_opt(:,k) = sol{2*na+1};
    rbar_old = rbar_opt(:,k);
  
    for i = 1 : na
        x_sol{i} = sol{i};
        u_sol{i} = sol{na+i};
    end

    % Opt
    for i = 1:na
        u_opt{i}(:,k) = u_sol{i}(:,1);
    end

    % Simulacao e evolucao dos estados
    for i = 1:na
    x{i}(:,k+1) =A*x{i}(:,k)+B*u_opt{i}(:,k);
    end
end

% tempo
t = 0.1*[0:kend];

%% Plots
figure (1)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),u_opt{1}(1,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t(1:end-1),u_opt{2}(1,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{3}(1,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{4}(1,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{5}(1,:),'r','LineStyle','--','LineWidth',2);
hold on; box on;
xlabel('Tempo [s]'); ylabel('u(1)')

figure (2)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),u_opt{1}(2,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t(1:end-1),u_opt{2}(2,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{3}(2,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{4}(2,:),'r','LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt{5}(2,:),'r','LineStyle','--','LineWidth',2);
hold on; box on;
xlabel('Tempo [s]'); ylabel('u(2)')

figure (3)
hold on
grid on; xlim([0 t(end)])
plot(t,x{1}(1,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t,x{2}(1,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{3}(1,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{4}(1,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{5}(1,:),'r','LineStyle','--','LineWidth',2)
hold on;
hold on; box on;
stairs(t(1:end-1),rbar_opt(1,:),'--k','linewidth',1)
xlabel('Tempo [s]')
ylabel('x1')
% legend('$x_1$','$x_2$','$\bar{R}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')

figure (4)
hold on
grid on; 
plot(t,x{1}(2,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t,x{2}(2,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{3}(2,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{4}(2,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{5}(2,:),'r','LineStyle','--','LineWidth',2)
hold on; box on;
stairs(t(1:end-1),rbar_opt(2,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo [s]')
ylabel('x2')
% legend('x_1','x_2')

figure (5)
hold on
grid on; 
plot(t,x{1}(3,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t,x{2}(3,:),'r','LineStyle',':','LineWidth',2)
plot(t,x{3}(3,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{4}(3,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{5}(3,:),'r','LineStyle','--','LineWidth',2)
hold on; box on;
hold on;xlim([0 t(end)])
xlabel('Tempo [s]')
ylabel('x3')

figure (6)
hold on
grid on; 
plot(t,x{1}(4,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t,x{2}(4,:),'r','LineStyle',':','LineWidth',2)
plot(t,x{3}(4,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{4}(4,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{5}(4,:),'r','LineStyle','--','LineWidth',2)
hold on; box on;
hold on;xlim([0 t(end)])
xlabel('Tempo [s]')
ylabel('x4')

figure (7)
hold on
grid on; 
plot(t,x{1}(5,:),'r','LineStyle','--','LineWidth',2); hold on;
plot(t,x{2}(5,:),'r','LineStyle',':','LineWidth',2)
plot(t,x{3}(5,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{4}(5,:),'r','LineStyle','--','LineWidth',2)
plot(t,x{5}(5,:),'r','LineStyle','--','LineWidth',2)
hold on; box on;
hold on;xlim([0 t(end)])
xlabel('Tempo [s]')
ylabel('x5')


eo =0;
for i = 1:5
    for j = 1:5
        if i==j
        else
            e = abs(x{i}(1,:)-x{j}(1,:))+abs(x{i}(2,:)-x{j}(2,:))+abs(x{i}(3,:)-x{j}(3,:))+abs(x{i}(4,:)-x{j}(4,:))+abs(x{i}(5,:)-x{j}(5,:))+eo;
            eo = e;
            
        end
    end
end
% c1 = sum(u_opt{});

figure (13)
hold on
grid on; xlim([0 t(end)])
plot(t,e,'--r','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('Somatorio dos erros entre os estados dos agentes')

% figure (14)
% hold on
% grid on; xlim([0 t(end)])
% plot(t(1:end-1),c,'--b','LineWidth',1.5); hold on;
% hold on;  box on;
% xlabel('Tempo [s]'); ylabel('Somatorio dos controles dos agentes')