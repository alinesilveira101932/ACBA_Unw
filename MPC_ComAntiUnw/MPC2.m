%% MPC - Usando YALMIP - Artigo MPC + Unwinding - Agentes separados
% Autor: Aline Isabel
% Data: 31/03/2026
clc;clear;
%close all;
%--------------------------------------------------------------------
% Load parametros
Parametros_yalmip
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45

yalmip('clear');
kend = 100;

x1(:,1) = [pi/3;10];
x2(:,1) = [pi/10;0];
x3(:,1) = [pi/3;-5];

for k = 1:kend

constraints = [];
objective = 0;

% 1. Inicializacao
Nbar1 = 0;
rbar = sdpvar(1,1); 
Nbar2 = intvar(1,1); Nbar3 = intvar(1,1);
U1 = sdpvar(1,N); U2 = sdpvar(1,N); U3 = sdpvar(1,N);
X1 = sdpvar(2,N+1);X2 = sdpvar(2,N+1); X3 = sdpvar(2,N+1);

% Equilibrio
RBAR1 = rbar + 2*pi*Nbar1; RBAR2 = rbar + 2*pi*Nbar2; RBAR3 = rbar + 2*pi*Nbar3;
xbar1 = Nx*RBAR1; xbar2 = Nx*RBAR2; xbar3 = Nx*RBAR3; 
ubar1 = Nu*RBAR1; ubar2 = Nu*RBAR2; ubar3 = Nu*RBAR3;

bf1 = bo - Or*RBAR1; bf2 = bo - Or*RBAR2; bf3 = bo - Or*RBAR3;

constraints = [constraints,X1(:,1) == x1(:,k)];
constraints = [constraints,X2(:,1) == x2(:,k)];
constraints = [constraints,X3(:,1) == x3(:,k)];

for i = 1:N
objective = objective + (X1(:,i) - xbar1)'*Q*(X1(:,i) - xbar1) + (U1(:,i) - ubar1)'*R*(U1(:,i) - ubar1);
objective = objective + (X2(:,i) - xbar2)'*Q*(X2(:,i) - xbar2) + (U2(:,i) - ubar2)'*R*(U2(:,i) - ubar2);
objective = objective + (X3(:,i) - xbar3)'*Q*(X3(:,i) - xbar3) + (U3(:,i) - ubar3)'*R*(U3(:,i) - ubar3);

constraints = [constraints,-5 <= U1(:,i) <= 5]; constraints = [constraints,-5 <= U2(:,i) <= 5]; constraints = [constraints,-5 <= U3(:,i) <= 5];
constraints = [constraints,X1(:,i+1)== A*X1(:,i)+B*U1(:,i)];
constraints = [constraints,X2(:,i+1)== A*X2(:,i)+B*U2(:,i)];
constraints = [constraints,X3(:,i+1)== A*X3(:,i)+B*U3(:,i)];

end
objective = objective + (X1(:,N+1) - xbar1)'*P*(X1(:,N+1) - xbar1);
objective = objective + (X2(:,N+1) - xbar2)'*P*(X2(:,N+1) - xbar2);
objective = objective + (X3(:,N+1) - xbar3)'*P*(X3(:,N+1) - xbar3);

constraints = [constraints, Sf*X1(:,N+1) <= bf1]; constraints = [constraints, Sf*X2(:,N+1) <= bf2]; constraints = [constraints, Sf*X3(:,N+1) <= bf3];

ops = sdpsettings('Verbose',0,'solver','gurobi');
solucao = optimize(constraints,objective,ops)

% Extrair valores ótimos
    x_opt1(:,k) = value(X1(:,1)); x_opt2(:,k) = value(X2(:,1)); x_opt3(:,k) = value(X3(:,1)); 
    u_opt1(:,k) = value(U1(:,1)); u_opt2(:,k) = value(U2(:,1)); u_opt3(:,k) = value(U3(:,1)); 
    rbar_opt(:,k) = value(rbar); N_opt2 = value(Nbar2); N_opt3 = value(Nbar3); 
    RBAR_opt(:,k) = rbar_opt(:,k) + 2*pi*[0; N_opt2; N_opt3];


    uk = [u_opt1(:,k); u_opt2(:,k); u_opt3(:,k)];
    xk = [x_opt1(:,k); x_opt2(:,k); x_opt3(:,k)];

   
xini1 = x1(:,k);
[t, xd1] = ode45(@(t,x) simulation1(x, u_opt1(:,k),Ac,Bc),[0 T], xini1, options2); 
x1(:,k+1) = xd1(length(t),:)';

xini2 = x2(:,k);
[t, xd2] = ode45(@(t,x) simulation2(x, u_opt2(:,k),Ac,Bc),[0 T], xini2, options2); 
x2(:,k+1) = xd2(length(t),:)';

xini3 = x3(:,k);
[t, xd3] = ode45(@(t,x) simulation3(x, u_opt3(:,k),Ac,Bc),[0 T], xini3, options2); 
x3(:,k+1) = xd3(length(t),:)';

end

function xdot = simulation1(x,u,Ac,Bc)
xdot = Ac*x+Bc*u;
end
function xdot = simulation2(x,u,Ac,Bc)
xdot = Ac*x+Bc*u;
end
function xdot = simulation3(x,u,Ac,Bc)
xdot = Ac*x+Bc*u;
end

t = 0.1*[0:kend];

%% Plots
figure (102)
hold on
grid on; xlim([0 t(end)])
plot(t(2:end),u_opt1,'LineStyle','--','Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t(2:end),u_opt2,'LineStyle','--','Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t(2:end),u_opt3,'LineStyle','--','Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u')

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x1(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x1(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
hold on; box on;
stairs(t(2:end),RBAR_opt(1,:),'--k','linewidth',1)
xlabel('Tempo, s')
ylabel('Agente 1')
legend('$x_1$','$x_2$','$\bar{R}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')

figure (7)
hold on
grid on; 
plot(t,x2(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x2(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),RBAR_opt(2,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 2')
% legend('x_1','x_2')

figure (8)
hold on
grid on; 
plot(t,x3(1,:),'LineStyle','--','Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x3(2,:),'LineStyle','--','Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t(2:end),RBAR_opt(3,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 3')


e12 = abs(exp(1j*x1(1,:))-exp(1j*x2(1,:))); e13 = abs(exp(1j*x1(1,:))-exp(1j*x3(1,:))); e23 = abs(exp(1j*x2(1,:))-exp(1j*x3(1,:)));

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