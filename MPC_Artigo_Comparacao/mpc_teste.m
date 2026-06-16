
clc;clear;close all;
%--------------------------------------------------------------------
% Load parametros
parametros_teste
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45

kend = 150;
 
% Estados iniciais

% x(:,1) = (-1+2*rand(25,1)); % numeros aleatorios no intervalo -1 a 1

controller = mpc_controller_teste(A,B,N,Q,R,S,na,nx,nu,I_na,beta); % Cria o problema

for k = 1:kend
    input = x(:,k);

    % Solucao
    [sol, diagnostic] = controller(input);
    u_sol = sol(N+2:end);
    x_sol = sol(1:N+1);
   
    % Opt
     u_opt(:,k) = u_sol{1};

    % Simulacao e evolucao dos estados
     x(:,k+1) = kron(I_na,A)*x(:,k)+kron(I_na,B)*u_opt(:,k);
   
end

% tempo
t = 0.1*[0:kend];

%% Plots
figure (102)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),u_opt(1,:),'LineStyle','--','LineWidth',2); hold on;
plot(t(1:end-1),u_opt(2,:),'LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt(3,:),'LineStyle','--','LineWidth',2)
plot(t(1:end-1),u_opt(4,:),'LineStyle','--','LineWidth',2)
plot(t(1:end-1),u_opt(5,:),'LineStyle','--','LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u1')
legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')

figure (103)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),u_opt(6,:),'LineStyle','--','LineWidth',2); hold on;
plot(t(1:end-1),u_opt(7,:),'LineStyle','--','LineWidth',2);
plot(t(1:end-1),u_opt(8,:),'LineStyle','--','LineWidth',2)
plot(t(1:end-1),u_opt(9,:),'LineStyle','--','LineWidth',2)
plot(t(1:end-1),u_opt(10,:),'LineStyle','--','LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u2')
legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x(1,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(2,:),'LineStyle','--','LineWidth',2)
plot(t,x(3,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(4,:),'LineStyle','--','LineWidth',2)
plot(t,x(5,:),'LineStyle','--','LineWidth',2); hold on;
hold on;
hold on; box on;
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')
xlabel('Tempo, s')
ylabel('x1')

figure (7)
hold on
grid on; xlim([0 t(end)])
plot(t,x(6,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(7,:),'LineStyle','--','LineWidth',2)
plot(t,x(8,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(9,:),'LineStyle','--','LineWidth',2)
plot(t,x(10,:),'LineStyle','--','LineWidth',2); hold on;
hold on;
hold on; box on;
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')
xlabel('Tempo, s')
ylabel('x2')

figure (8)
hold on
grid on; xlim([0 t(end)])
plot(t,x(11,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(12,:),'LineStyle','--','LineWidth',2)
plot(t,x(13,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(14,:),'LineStyle','--','LineWidth',2)
plot(t,x(15,:),'LineStyle','--','LineWidth',2); hold on;
hold on;
hold on; box on;
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')
xlabel('Tempo, s')
ylabel('x3')

figure (9)
hold on
grid on; xlim([0 t(end)])
plot(t,x(16,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(17,:),'LineStyle','--','LineWidth',2)
plot(t,x(18,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(19,:),'LineStyle','--','LineWidth',2)
plot(t,x(20,:),'LineStyle','--','LineWidth',2); hold on;
hold on;
hold on; box on;
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')
xlabel('Tempo, s')
ylabel('x4')

figure (10)
hold on
grid on; xlim([0 t(end)])
plot(t,x(21,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(22,:),'LineStyle','--','LineWidth',2)
plot(t,x(23,:),'LineStyle','--','LineWidth',2); hold on;
plot(t,x(24,:),'LineStyle','--','LineWidth',2)
plot(t,x(25,:),'LineStyle','--','LineWidth',2); hold on;
hold on;
hold on; box on;
hold on; box on; legend('Ag1', 'Ag2', 'Ag3','Ag4','Ag5')
xlabel('Tempo, s')
ylabel('x5')

eo =0;
for i = 1:5
    for j = 1:5
        if i==j
        else
            e = abs(x(i,:)-x(j,:))+abs(x(i+5,:)-x(j+5,:))+abs(x(i+10,:)-x(j+10,:))+abs(x(i+15,:)-x(j+15,:))+abs(x(i+20,:)-x(j+20,:))+eo;
            eo = e;
            
        end
    end
end
c = sum(u_opt);

figure (13)
hold on
grid on; xlim([0 t(end)])
plot(t,e,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('Somatorio dos erros entre os estados dos agentes')

figure (14)
hold on
grid on; xlim([0 t(end)])
plot(t(1:end-1),c,'--b','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo [s]'); ylabel('Somatorio dos controles dos agentes')