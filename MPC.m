%% Controle Preditivo - Artigo  MPC + Unwinding
% Codigo geral
% Aline
% 23/01/2026
%-----------------------------------
clc; clear; close all;

% Load Parametros da planta
Parametros

%% Determinar O_psi_inf - todas as refs

Gamma = [eye(nx*na) zeros(nx*na,nr); -K (K*Nx+Nu); zeros(na*nx,nx*na) Nx; zeros(nu*na,nx*na) Nu];

Spsi = blkdiag(Sx_barra,Su_barra,Sx_barra,Su_barra);
bpsi = [bx_barra; bu_barra; bx_barra-ex; bu_barra-eu];

Apsi_f = [Af B*(K*Nx+Nu); zeros(nq*nr,nx*na) eye(nq*nr)];
max_iter = 50;

[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

% figure (1)
% plot(O_psi_inf)

%% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx*na);
Or = So(:,nx*na+1:end);
Sf = Ox;

%% Dual Mode Predictive Control

% Carregando matrizes de parametros do controle preditivo
Matrizes_Predicao

% Configurando solvers
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45
options_qudprog = optimoptions('quadprog','Display','off');

% Inicializacao
kend = 100;

% Loop
x(:,1) = [pi/3;10;pi/10;0;pi/3;-5];
x0 = x(:,1);
sys = ss(Ac,Bc,[],[]);
trec = 0; xrec = x(:,1)'; nt = 10; 

rbar_opt(:,1) = [0;0;0];


    u_var = sdpvar(N*nu*na,1);
    rbar_var = sdpvar(1,1);
    N_var = intvar(2,1);

for k = 1:kend

  % Matrizes finais do quadprog 
    % J = f*U'+0.5*U'*H*U
    % s.a A * U < B
     % 1. Inicializacao

    rbar = rbar_var*ones(nr,1) + 2*pi*[0; N_var];

    dummy = [u_var; rbar];

    bf = bo - Or*rbar;
    
    bxp = [repmat(bx_barra,N-1,1);bf];

    Aqp = [Sxp*Bp   zeros(size(Sxp*Bp,1),nr);  Sup zeros(size(Sup,1),nr)];
    bqp = [bxp - Sxp*Ap*x(:,k);       bup];

    Hqp = [Bp'*Qp*Bp + Rp,              (-(Inx*Nx)'*Qp*Bp-(Inu*Nu)'*Rp)';
          -((Inx*Nx)'*Qp*Bp-(Inu*Nu)'*Rp),      (Nx'*(Q_barra + Inx'*Qp*Inx)*Nx + Nu'*(Inu'*Rp*Inu)*Nu)]; 
        
    
    % Simetrização
    Hqp = (Hqp + Hqp')/2;

    fqp = [(Ap*x(:,k))'*Qp*Bp,           ((-x(:,k)'*Q_barra-(Ap*x(:,k))'*Qp*Inx)*Nx)] ;

    
    obj = 0.5*dummy'*Hqp*dummy+fqp*dummy;

    LMIs = [];

    % 3. LMIs
    LMIs = [LMIs,Aqp*dummy-bqp <= 0];

    % 4. Resolver as LMIs 
    solucao = optimize(LMIs, obj, sdpsettings('solver','gurobi'));
    if solucao.problem ~= 0
        solucao.info
    end
    
    dummy_opt = value(dummy);
    u(:,k) = dummy_opt(1:na);
  
    rbar_1(:,k) = dummy_opt(na*N*nu+1);
    rbar_2(:,k) = dummy_opt(na*N*nu+2);
    rbar_3(:,k) = dummy_opt(end);
    rbar_opt(:,k+1) = [rbar_1(:,k);rbar_2(:,k);rbar_3(:,k)];

    N_var_opt(:,k) = value(N_var);
    rbar_var_opt(:,k) = value(rbar_var);

    xini = x(:,k);
    [t, xd] = ode45(@(t,x) simulation(x, u(:,k),Ac_barra,Bc_barra),[0 T], xini, options2); 
    x(:,k+1) = xd(length(t),:)';
    
    % 
    % t = linspace((k-1)*T,k*T,nt);
    % U = repmat(u(:,k)',nt,1);
    % [~,tout,xout] = lsim(sys,U,t,x(:,k));
    % trec = [trec;tout(2:end)];
    % xrec = [xrec;xout(2:end,:)];
    % x(:,k+1) = xrec(end,:)';

end

function xdot = simulation(x,u,Ac,Bc)
xdot = Ac*x+Bc*u;
end

%% Plots
t = 0.1*(0:kend);

figure (6)
hold on
grid on; xlim([0 t(end)])
plot(t,x(1,:),'Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(2,:),'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on;
hold on; box on;
stairs(t,rbar_opt(1,:),'--k','linewidth',1)
xlabel('Tempo, s')
ylabel('Agente 1')
legend('$x_1$','$x_2$','$\bar{R}_i$', 'Interpreter', 'latex', 'FontSize', 14,'Location','northoutside','Orientation','horizontal')

figure (7)
hold on
grid on; 
plot(t,x(3,:),'Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(4,:),'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t,rbar_opt(2,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 2')
% legend('x_1','x_2')

figure (8)
hold on
grid on; 
plot(t,x(5,:),'Color',[0.85,0.33,0.10],'LineWidth',2); hold on;
plot(t,x(6,:),'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on; box on;
stairs(t,rbar_opt(3,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 t(end)])
xlabel('Tempo, s')
ylabel('Agente 3')


figure (2)
hold on
grid on; xlim([0 t(end)])
plot(t(2:end),u(1,:),'Color',[0.39,0.83,0.07],'LineWidth',2); hold on;
plot(t(2:end),u(2,:),'Color',[0.07,0.62,1.00],'LineWidth',2);
plot(t(2:end),u(3,:),'Color',[0.93,0.69,0.13],'LineWidth',2)
hold on; box on;
xlabel('Tempo, s'); ylabel('u')
legend('Agente 1','Agente 2', 'Agente 3')
plot([0 t(end)],[umax1 umax1],'--r','linewidth',1,'handlevisibility','off')
plot([0 t(end)],[umin1 umin1],'--r','linewidth',1,'handlevisibility','off')
text(0.1,5.3,'Limitante de controle','Color','r','FontWeight','bold')
ylim([umin1-0.5 umax1+0.5])


figure (104)
subplot(4,1,1)
grid on; 
plot([0 t(end)],[0 0],'Color',[0.39,0.83,0.07],'LineWidth',2)
xlim([0 t(end)])
xlabel('Tempo, s');  box on;
legend('N_1')

subplot(4,1,2)
stairs(t(2:end),N_var_opt(1,:),'Color',[0.07,0.62,1.00],'LineWidth',2); hold on;
xlim([0 t(end)]); grid on; box on;
xlabel('Tempo, s');  
legend('N_2')

subplot(4,1,3)
stairs(t(2:end),N_var_opt(2,:),'Color',[0.93,0.69,0.13],'LineWidth',2); hold on;grid on;
xlim([0 t(end)])
xlabel('Tempo, s');   box on;
legend('N_3')

subplot(4,1,4)
plot(t(2:end),rbar_var_opt,'k','LineWidth',2); hold on;grid on;
xlim([0 t(end)])
xlabel('Tempo, s');   box on;
legend('$\bar{r}$', 'Interpreter', 'latex')


%% Resultados

D2pi_com = (mod(x(1,:)-x(3,:)+pi,2*pi)-pi).^2 + abs(rem(x(1,:)-x(5,:),2*pi)).^2 + abs(rem(x(5,:)-x(3,:),2*pi)).^2 + abs(x(2,:)-x(4,:)).^2+abs(x(2,:)-x(6,:)).^2+abs(x(6,:)-x(4,:)).^2;
check1 = abs(rbar_opt(1,end)- rbar_opt(2,end))/(2*pi)
check2 = abs(rbar_opt(1,end)- rbar_opt(3,end))/(2*pi)
check3 = abs(rbar_opt(3,end)- rbar_opt(2,end))/(2*pi)

% figure (11)
% hold on
% grid on; xlim([0 t(end)])
% plot(t,D2pi_com,'r','LineWidth',2); hold on;
% hold on;  box on;
% xlabel('Tempo, s'); ylabel('D')
% legend('Com Anti-Unwinding','Sem Anti-Unwinding')

e12 = abs(exp(1j*x(1,:))-exp(1j*x(3,:))); e13 = abs(exp(1j*x(1,:))-exp(1j*x(5,:))); e23 = abs(exp(1j*x(3,:))-exp(1j*x(5,:)));


figure (13)
subplot(3,1,1)
hold on
grid on; xlim([0 t(end)])
plot(t,e12,'-r','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{12}')

subplot(3,1,2)
hold on
grid on; xlim([0 t(end)])
plot(t,e13,'--r','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{13}')
subplot(3,1,3)
hold on
grid on; xlim([0 t(end)])
plot(t,e23,'-.r','LineWidth',1.5); hold on;
hold on;  box on;
xlabel('Tempo, s'); ylabel('e_{23}')