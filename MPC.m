%% Controle Preditivo - Artigo  MPC + Unwinding
% Codigo geral
% Aline
% 23/01/2026
%-----------------------------------
clc; clear; close all;

% Load Parametros da planta
Parametros

%% Determinar O_psi_inf - todas as refs

Gamma = [eye(nx*na) zeros(nx*na,1); -K (K*Nx+Nu); zeros(na*nx,nx*na) Nx; zeros(nu*na,nx*na) Nu];

Spsi = blkdiag(Sx_barra,Su_barra,Sx_barra,Su_barra);
bpsi = [bx_barra; bu_barra; bx_barra-ex; bu_barra-eu];

Apsi_f = [Af B*(K*Nx+Nu); zeros(nq,nx*na) eye(nq)];
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
x(:,1) = [5;0;0.2;0;0.15;0];

sys = ss(Ac,Bc,[],[]);
trec = 0; xrec = x(:,1)'; nt = 10; 

rbar(:,1) = 0;


for k = 1:kend

  % Matrizes finais do quadprog 
    % J = f*U'+0.5*U'*H*U
    % s.a A * U < B

    bf = bo - Or*rbar(:,k);
    
    bxp = [repmat(bx_barra,N-1,1);bf];

    Aqp = [Sxp*Bp   zeros(615,1);  Sup zeros(150,1)];
    bqp = [bxp - Sxp*Ap*x(:,k);       bup];

    Hqp = [Bp'*Qp*Bp + Rp,              (-(Inx*Nx)'*Qp*Bp-(Inu*Nu)'*Rp)';
          -((Inx*Nx)'*Qp*Bp-(Inu*Nu)'*Rp),      (Nx'*(Q_barra + Inx'*Qp*Inx)*Nx + Nu'*(Inu'*Rp*Inu)*Nu)]; 
        
    
    % Simetrização
    Hqp = (Hqp + Hqp')/2;

    fqp = [(Ap*x(:,k))'*Qp*Bp,           ((-x(:,k)'*Q_barra-(Ap*x(:,k))'*Qp*Inx)*Nx)] ;

    dummy = quadprog(Hqp, fqp, [],[], [], [], [], [], [], options_qudprog);

    u(:,k) = dummy(1:3);
    rbar_opt(:,k) = dummy(end);


    % simulacao linear
    % x(:,k+1) = A*x(:,k)+B*u(:,k);
    
    rbar(:,k+1) = rbar_opt(:,k);
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
% w = sqrt((a*sin(x(1))+u)/c);
% xdot_1 = x(2);
% xdot_2 = -a*sin(x(1))-b*x(2)+c*1.05*w^2;
% xdot = [xdot_1;xdot_2];
xdot = Ac*x+Bc*u;
end

% Plots
figure (1)
hold on
grid on; 
plot(u(1,:),'r','LineWidth',2)
hold on;
title('Agente1')
xlabel('Time, s')
ylabel('u, rad')

figure (2)
hold on
grid on; 
plot(u(2,:),'r','LineWidth',2)
hold on;
title('Agente2')
xlabel('Time, s')
ylabel('u, rad')

figure (4)
hold on
grid on; 
plot(u(3,:),'r','LineWidth',2)
hold on;
title('Agente3')
xlabel('Time, s')
ylabel('u, rad')

figure (6)
hold on
grid on; 
plot(x(1,:),'r','LineWidth',2); hold on;
plot(x(2,:),'b','LineWidth',2)
hold on;
hold on;
xlabel('k')
title('Agente1')
legend('x1','x2')

figure (7)
hold on
grid on; 
plot(x(3,:),'r','LineWidth',2); hold on;
plot(x(4,:),'b','LineWidth',2)
hold on;
hold on;
xlabel('k')
title('Agente2')
legend('x1','x2')

figure (8)
hold on
grid on; 
plot(x(5,:),'r','LineWidth',2); hold on;
plot(x(6,:),'b','LineWidth',2)
hold on;
hold on;
xlabel('k')
title('Agente3')
legend('x1','x2')

figure (9)
hold on
grid on; 
plot(rbar,'r','LineWidth',2); hold on;
hold on;
hold on;
xlabel('k')
legend('rbar')