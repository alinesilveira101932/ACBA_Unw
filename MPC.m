%% MPC - Usando YALMIP - Artigo MPC + Unwinding - Agentes separados
% Autor: Aline Isabel
% Data: 31/03/2026
clc;clear;
%close all;
%--------------------------------------------------------------------
% Load parametros
Parametros_yalmip
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45


kend = 100;
 
% Estados iniciais
x{1}(:,1) = [pi/3;10];
x{2}(:,1) = [pi/10;0];
x{3}(:,1) = [pi/3;-5];

controller = mpc_controller(A,B,N,Q,R,P,Nx,Nu,bo,Or,Sf,na,nx,nu); % Cria o problema

for k = 1:kend
    input = [x{1}(:,k); x{2}(:,k); x{3}(:,k)];

    % Solucao
    sol = controller(input);

    rbar_opt(:,k) = sol{2*na+1};
    Nbar_opt{1}(:,k) = 0;
    RBAR_opt{1}(:,k) = rbar_opt(:,k);
    for i = 1: na-1
        Nbar_opt{i+1}(:,k) = sol{2*na+1+i};
        RBAR_opt{i+1}(:,k) = rbar_opt(:,k) + 2*pi*Nbar_opt{i+1}(:,k);
    end
    for i = 1 : na
        x_sol{i} = sol{i};
        u_sol{i} = sol{na+i};
    end

    % Opt
    for i = 1:na
        u_opt{i}(:,k) = u_sol{i}(1);
    end

    % Simulacao e evolucao dos estados
    for i = 1:na
    xini = x{i}(:,k);
    [t, xd] = ode45(@(t,x) simulation(x, u_opt{i}(:,k),Ac,Bc),[0 T], xini, options2); 
    x{i}(:,k+1) = xd(length(t),:)';
    end
end

function xdot = simulation(x,u,Ac,Bc)
xdot = Ac*x+Bc*u;
end


% tempo
t = 0.1*[0:kend];


%% Plots
figure (102)
hold on
grid on; xlim([0 t(end)])
plot(t(2:end),u_opt{1},'LineStyle','--','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t(2:end),u_opt{2},'LineStyle','--','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t(2:end),u_opt{3},'LineStyle','--','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u')

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x{1}(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x{1}(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
hold on; box on;
stairs(t(2:end),RBAR_opt{1},'--k','linewidth',1)
xlabel('Tempo, s')
ylabel('Agente 1')
legend('$x_1$','$x_2$','$\bar{R}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')

figure (7)
hold on
grid on; 
plot(t,x{2}(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x{2}(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),RBAR_opt{2},'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 2')
% legend('x_1','x_2')

figure (8)
hold on
grid on; 
plot(t,x{3}(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x{3}(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),RBAR_opt{3},'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 3')


e12 = abs(exp(1j*x{1}(1,:))-exp(1j*x{2}(1,:))); e13 = abs(exp(1j*x{1}(1,:))-exp(1j*x{3}(1,:))); e23 = abs(exp(1j*x{2}(1,:))-exp(1j*x{3}(1,:)));

figure (13)
subplot(3,1,1)
hold on
grid on; xlim([0 t(end)])
plot(t,e12,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{12}')

subplot(3,1,2)
hold on
grid on; xlim([0 t(end)])
plot(t,e13,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{13}')
subplot(3,1,3)
hold on
grid on; xlim([0 t(end)])
plot(t,e23,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{23}')