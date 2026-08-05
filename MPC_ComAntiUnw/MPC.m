%% MPC - Usando YALMIP - COM anti Unwinding - Agentes separados
% Custo novo: reformulado para usar xbar^* e ubar^*, nao mais rbar*
% codigo generalizado para Na agentes
% agora usa uma parcela ||xbar(k)-xbar(k-1)||_W
% Autor: Aline Isabel
% Data: 05/08/2026
clc;clear;
%close all;
%--------------------------------------------------------------------
% Load parametros
Parametros_yalmip
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45

tic
Time =  tic;

kend = 100;
 
% Estados iniciais - Int Duplo - ocorre unw
x{1}(:,1) = [pi/3;10]; x{2}(:,1) = [pi/10;0]; x{3}(:,1) = [pi/3;-5]; x{4}(:,1) = [pi/4;-1];
x{5}(:,1) = [pi/(3.2);8];
x{6}(:,1) = [pi/2;15];
x{7}(:,1) = [0;0];
x{8}(:,1) = [3*pi/2;-4];
x{9}(:,1) = [3*pi/2;8];
x{10}(:,1) = [2*pi-0.3;3];
x{11}(:,1) = [pi/2;-15];
x{12}(:,1) = [pi;0];
x{13}(:,1) = [pi;20];
x{14}(:,1) = [0;-17];
x{15}(:,1) = [pi/5;3];
x{16}(:,1) = [pi/(7);8];
x{17}(:,1) = [pi/20;10];
x{18}(:,1) = [3*pi/5;10];
x{19}(:,1) = [3*pi/2;1];
x{20}(:,1) = [pi-0.3;1];

% Estados iniciais - Int Duplo - nao ocorre unw
% x{1}(:,1) = [0;0];
% x{2}(:,1) = [pi/5;0];
% x{3}(:,1) = [pi/2;0];

for i = 1 :na
XBAR_old{i} = zeros(nx,1);
end

controller = mpc_controller(A,B,N,Q,R,P,W,bo,Oxbar,Oubar,Sf,na,nx,nu,Su,bu); % Cria o problema

for k = 1:kend

    % inicial
    input = [];
    for i = 1:na
        input = [input ;x{i}(:,k)];
    end
     for i = 1:na
        input = [input ;XBAR_old{i}];
    end
   
    % Solucao
    [sol, diagnostics] = controller(input);
    Nbar_opt{1}(:,k) = [0;0];
    
    for i = 1: na-1
        Nbar_opt{i+1}(:,k) = sol{4*na+i};
    end
    for i = 1 : na
        x_sol{i} = sol{i};
        u_sol{i} = sol{na+i};
        XBAR_opt{i}(:,k) = sol{2*na+i};
        UBAR_opt{i}(:,k) = sol{3*na+i};
        XBAR_old{i} =  XBAR_opt{i}(:,k);
    end
    
    % Opt
    for i = 1:na
        u_opt{i}(:,k) = u_sol{i}(1);
    end

% Simulacao e evolucao dos estados
    for i = 1:na
    x{i}(:,k+1) =A*x{i}(:,k)+B*u_opt{i}(:,k);
    end

    % for i = 1:na
    % custo = custo + (x{i}(:,k)-xbar_opt(:,k))'*Q*(x{i}(:,k)-xbar_opt(:,k))+(u_opt{i}(:,k)-ubar_opt(:,k))'*R*(u_opt{i}(:,k)-ubar_opt(:,k));
    % end
end

% tempo
toc
Tempo_total = toc(Time)
t = 0.1*[0:kend];


%% Plots
figure (102)
hold on
grid on; xlim([0 t(end)])
for i = 1:na
plot(t(1:end-1),u_opt{i},'LineStyle','-','LineWidth',2); hold on;
end
hold on; box on;
xlabel('Tempo [s]'); ylabel('u')


for i = 1:na
figure (i)
hold on
grid on; xlim([0 t(end)])
plot(t,x{i}(1,:),'LineStyle','-','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x{i}(2,:),'LineStyle','-','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
stairs(t(1:end-1),XBAR_opt{i}(1,:),'--k','linewidth',1)
hold on; box on;
xlabel('Tempo [s]')
ylabel(sprintf('Agente %d', i))
legend('$x_1$','$x_2$','$\bar{x}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')
end

figure (100)
hold on
grid on; xlim([0 t(end)])
for i = 1:na
plot(t,x{i}(1,:),'LineStyle','-','LineWidth',2); hold on;
end
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5','Ag6','Ag7','Ag8','Ag9','Ag10','Ag11','Ag12','Ag13','Ag14','Ag15','Ag16','Ag17','Ag18','Ag19','Ag20')
xlabel('Tempo [s]'); ylabel('x1')

figure (101)
hold on
grid on; xlim([0 t(end)])
for i =1:na
plot(t,x{i}(2,:),'LineStyle','-','LineWidth',2); hold on;
end
hold on; box on; legend('Ag1', 'Ag2', 'Ag3')
xlabel('Tempo [s]'); ylabel('x2')


e12 = abs(exp(1j*x{1}(1,:))-exp(1j*x{2}(1,:))); e13 = abs(exp(1j*x{1}(1,:))-exp(1j*x{3}(1,:))); e23 = abs(exp(1j*x{2}(1,:))-exp(1j*x{3}(1,:)));

figure (13)
subplot(3,1,1)
hold on
grid on; xlim([0 t(end)])
plot(t,e12,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('e_{12}')

subplot(3,1,2)
hold on
grid on; xlim([0 t(end)])
plot(t,e13,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('e_{13}')
subplot(3,1,3)
hold on
grid on; xlim([0 t(end)])
plot(t,e23,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('e_{23}')
% 
% f12 = abs(exp(1j*x{1}(2,:))-exp(1j*x{2}(2,:))); f13 = abs(exp(1j*x{1}(2,:))-exp(1j*x{3}(2,:))); f23 = abs(exp(1j*x{2}(2,:))-exp(1j*x{3}(2,:)));
% 
% figure (14)
% subplot(3,1,1)
% hold on
% grid on; xlim([0 t(end)])
% plot(t,f12,'--b','LineWidth',1.5); hold on;
% hold on;  box on;
% xlabel('Tempo [s]'); ylabel('f_{12}')
% 
% subplot(3,1,2)
% hold on
% grid on; xlim([0 t(end)])
% plot(t,f13,'--b','LineWidth',1.5); hold on;
% hold on;  box on;
% xlabel('Tempo [s]'); ylabel('f_{13}')
% subplot(3,1,3)
% hold on
% grid on; xlim([0 t(end)])
% plot(t,f23,'--b','LineWidth',1.5); hold on;
% hold on;  box on;
% xlabel('Tempo [s]'); ylabel('f_{23}')
% 
e = zeros(1,size(x{1},2));   % 1 x tempo

for i = 1:na-1
    for j = i+1:na
        e = e + sum(abs(exp(1j*x{i}(1,:))-exp(1j*x{j}(1,:))),1)+ sum(abs(exp(1j*x{i}(2,:))-exp(1j*x{j}(2,:))),1);
    end
end

% c1 = sum(u_opt{});

figure (130)
hold on
grid on; xlim([0 t(end)])
plot(t,e,'--r','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('Soma dos erros de estado')