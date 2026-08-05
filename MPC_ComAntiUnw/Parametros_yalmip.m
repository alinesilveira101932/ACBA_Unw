%% Parametros - Usando YALMIP - Artigo MPC + Unwinding - Agentes separados
% Auxiliar
% Autor: Aline Isabel
% Data: 31/03/2026

%--------------------------------------------------------------------
% Planta
% Ac = [0 1; 0 0];
% Bc = [0;1];
% C = [1 0];

b = 0.5;
Ac = [0 1; 0 -b];
Bc = [0;1];
C = [1 0];

T = 0.1; % Periodo de amostragem

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
nq = size(C,1);

na = 20; % numero de agentes

N = 30; % Horizonte de predicao


% Matrizes de Peso
% Q = [3 0;0 10];
% R = 1;
% W = 0.01;
Q = 10*eye(nx);
R = 1;
W = 0.1; 

% Discretizacao
[A,B] = c2dm(Ac,Bc,[],[],T, 'zoh');

% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;

% Valores de equilibrio
aux1 = inv([A-eye(nx) B;C zeros(1)])*[zeros(nx,1);eye(1)];
Nx = aux1(1:nx,:); Nu = aux1(nx+1:end,:);

%% Restricoes 

% Su *u <= bu
umax = 5;  umin = -5;

% Sx*x <= bx
x1max = 1000; x1min = -1000;
x2max = 1000; x2min = -1000;

Su = [1; -1];
eu = 1e-3*ones(1*size(Su,1),1);

Sx = [1 0; -1 0; 0 1; 0 -1]; 
ex = 1e-3*ones(1*size(Sx,1),1);

bx = [x1max ; -x1min; x2max; -x2min];
bu = [umax; -umin];

Sxbar = [A-eye(nx);A-eye(nx)];
bxubar = [zeros(nx,1); zeros(nx,1)];
Subar = [B;B];


% %% Determinar O_psi_inf - com respeito a xbar -> sendo ubar = 0;
% 
% Gamma = [eye(nx) zeros(nx); -K K; zeros(nx,nx) eye(nx)];
% Spsi = blkdiag(Sx,Su,[Sx ;Sxbar]);
% bpsi = [bx; bu; [bx-ex;bxbar]];
% 
% Apsi_f = [Af B*K; zeros(nx,nx) eye(nx)];
% max_iter = 100;
% 
% [So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter);
% 
% O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox
% 
% % figure (1)
% % plot(O_psi_inf)
% 
% % Determinar Oinf a partir de O_psi_inf
% Ox = So(:,1:nx);
% Oxbar = So(:,nx+1:end);
% Sf = Ox;

%% Determinar O_psi_inf - com respeito a xbar -> sendo ubar = 0;

Gamma = [eye(nx) zeros(nx) zeros(nx,nu); -K K eye(nu); zeros(nx,nx) eye(nx) zeros(nx,nu); zeros(nu,nx) zeros(nu,nx) eye(nu)];
Spsi_1 = blkdiag(Sx,Su,Sx,Su);
bpsi_1 = [bx; bu; bx-ex;bu-eu];

Spsi = [Spsi_1; [zeros(2*nx,nx+nu) Sxbar Subar]];
bpsi = [bx; bu; bx-ex;bu-eu;bxubar];
Apsi_f = [Af B*K B; zeros(nx,nx) eye(nx) zeros(nx,nu); zeros(nu,nx) zeros(nu,nx) eye(nu)];
max_iter = 100;
tol = 1e-4;
[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter,tol);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

% figure (1)
% plot(O_psi_inf)

% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Oxbar = So(:,nx+1:nx*nx);
Oubar = So(:,nx*nx+1:end);
Sf = Ox;
