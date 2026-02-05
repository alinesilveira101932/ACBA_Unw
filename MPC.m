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
x(:,1) = [10;0;0.3;0;-1;0];
x0 = x(:,1);
sys = ss(Ac,Bc,[],[]);
trec = 0; xrec = x(:,1)'; nt = 10; 

rbar_opt(:,1) = [0;0;0];
M(:,1) = [0;0;0];
rbar_var0 = 0; 
N_var0 = [0;0];
u_var0 = zeros(na*nu*N,1);

for k = 1:kend

  % Matrizes finais do quadprog 
    % J = f*U'+0.5*U'*H*U
    % s.a A * U < B

    % 1. Inicializacao
    u_var = sdpvar(N*nu*na,1);
    rbar_var = sdpvar(1,1);
    N_var = intvar(2,1);
    
    % Warm start
    assign(rbar_var, rbar_var0);
    assign(N_var, N_var0);
    % assign(u_var, u_var0)

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
    solucao = optimize(LMIs, obj, sdpsettings('solver','gurobi','usex0',1));
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
    
    % Warm start
    rbar_var0 = rbar_var_opt(:,k);
    N_var0 =  N_var_opt(:,k) ;
    rbar0 = rbar_var0*ones(nr,1) + 2*pi*[0; N_var0];

    UN_1 = dummy_opt(1:(N-1)*na);


    xN = Ap*x(:,k) + Bp*dummy_opt(1:end-nr);
    xN_1 = xN((N-2)*na*nx+1:(N-1)*na*nx);
 
    u_var0 = [UN_1; -K*(xN_1-Nx*rbar0)+Nu*rbar0];

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

% Plots
figure (1)
hold on
grid on; 
plot(u(1,:),'r','LineWidth',2)
hold on;
title('Agente1')
xlabel('k')
ylabel('u, rad')
hold on;
% plot([0 kend],[umax1 umax1],'--k','linewidth',1,'handlevisibility','off')
% text(0.1,0.18,'Limitante de x2','Color','r','FontWeight','bold')
xlim([0 kend])

figure (2)
hold on
grid on; 
plot(u(2,:),'r','LineWidth',2)
hold on; xlim([0 kend])
title('Agente2')
xlabel('k')
ylabel('u, rad')
% plot([0 kend],[umin2 umin2],'--k','linewidth',1,'handlevisibility','off')

figure (4)
hold on
grid on; 
plot(u(3,:),'r','LineWidth',2)
hold on;xlim([0 kend])
title('Agente3')
xlabel('k')
ylabel('u, rad')
% plot([0 kend],[umin3 umin3],'--k','linewidth',1,'handlevisibility','off')

figure (6)
hold on
grid on; xlim([0 kend])
plot(x(1,:),'r','LineWidth',2); hold on;
plot(x(2,:),'b','LineWidth',2)
hold on;
hold on;
stairs(rbar_opt(1,:),'--k','linewidth',1,'handlevisibility','off')
xlabel('k')
title('Agente1')
legend('x1','x2')

figure (7)
hold on
grid on; 
plot(x(3,:),'r','LineWidth',2); hold on;
plot(x(4,:),'b','LineWidth',2)
hold on;
stairs(rbar_opt(2,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 kend])
xlabel('k')
title('Agente2')
legend('x1','x2')

figure (8)
hold on
grid on; 
plot(x(5,:),'r','LineWidth',2); hold on;
plot(x(6,:),'b','LineWidth',2)
hold on;
stairs(rbar_opt(3,:),'--k','linewidth',1,'handlevisibility','off')
hold on;xlim([0 kend])
xlabel('k')
title('Agente3')
legend('x1','x2')

figure (9)
hold on
grid on; 
stairs(rbar_opt(1,:),'LineWidth',2); hold on;
stairs(rbar_opt(2,:),'LineWidth',2); hold on;
stairs(rbar_opt(3,:),'LineWidth',2); hold on;
hold on;
hold on;xlim([0 kend])
xlabel('k')
legend('rbar1','rbar2','rbar3')

figure (10)
hold on
grid on; 
plot([0 kend],[0 0],'LineWidth',2)
stairs(N_var_opt(1,:),'LineWidth',2); hold on;
stairs(N_var_opt(2,:),'LineWidth',2); hold on;

plot(rbar_var_opt,'LineWidth',2); hold on;
hold on;
hold on;xlim([0 kend])
xlabel('k')
legend('N1','N2','N3','rbar')

check1 = abs(rbar_opt(1,end)- rbar_opt(2,end))/(2*pi)
check2 = abs(rbar_opt(1,end)- rbar_opt(3,end))/(2*pi)
check3 = abs(rbar_opt(3,end)- rbar_opt(2,end))/(2*pi)
