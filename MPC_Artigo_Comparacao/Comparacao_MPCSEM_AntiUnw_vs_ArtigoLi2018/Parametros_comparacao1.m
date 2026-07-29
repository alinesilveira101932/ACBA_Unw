%% Parametros - codigo auxiliar
% MPC - Comparacao sem anti Unwinding  vs artigo li2018
% Modelo artigo Li2018 - 5 estados e 2 entradas
% Autor: Aline Isabel
% Data: 24/07/2026
%--------------------------------------------------------------------
A = [0.8,0.1,0.1,0,0;
    0,0.9,0,0.1,0;
    0.1,0.1,0.6,0.1,0.1;
    0,0.1,0.1,0.8,0;
    0.1,0.1,0,0,0.8];

B = [-0.1,0.1;
    0.1,-0.2;
    0,-0.3;
    0.08,0.1;
    0.2,0.08];
C = [1 0 0 0 0;0 1 0 0 0]; 

N = 15; 

% Dimencoes da planta
nx = length(A);
nu = size(B,2);
ny = size(C,1);

na = 5; % numero de agentes

%% Tecnica MPC sem anti-unw
% Matrizes de Peso
Q = eye(5); R = eye(2);
W = 0.01;

% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;

% Valores de equilibrio
aux1 = inv([A-eye(nx) B;C zeros(ny,nu)])*[zeros(nx,nu);eye(ny,nu)];
Nx = aux1(1:nx,:); Nu = aux1(nx+1:end,:);

%% Restricoes 
umax = 0.3;  umin = -0.3;
Su = [1 0; -1 0;0 1;0 -1];
bu = [umax; -umin;umax;-umin];
eu = 1e-3*ones(1*size(bu,1),1);

% Sx*x <= bx
xmax = 10000; xmin = -10000;


Sx = [1 0 0 0 0; -1 0 0 0 0; 0 1 0 0 0; 0 -1 0 0 0; 0 0 1 0 0; 0 0 -1 0 0; 0 0 0 1 0; 0 0 0 -1 0;0 0 0 0 1;0 0 0 0 -1]; 
ex = 1e-3*ones(1*size(Sx,1),1);

bx = [xmax ; -xmin; xmax; -xmin;xmax; -xmin;xmax; -xmin;xmax; -xmin];


%% Determinar O_psi_inf - todas as refs

Gamma = [eye(nx) zeros(nx,ny); -K (K*Nx+Nu); zeros(nx,nx) Nx; zeros(nu,nx) Nu];

Spsi = blkdiag(Sx,Su,Sx,Su);
bpsi = [bx; bu; bx-ex; bu-eu];

Apsi_f = [Af B*(K*Nx+Nu); zeros(ny,nx) eye(ny)];
max_iter = 100;

[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

% figure (1)
% plot(O_psi_inf)

% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Or = So(:,nx+1:end);
Sf = Ox;


%% Determinar O_psi_inf - todas as refs

% Gamma = [-K (K*Nx+Nu); zeros(nu,nx) Nu];
% 
% Spsi = blkdiag(Su,Su);
% bpsi = [bu; bu-eu];
% 
% Apsi_f = [Af B*(K*Nx+Nu); zeros(ny,nx) eye(ny)];
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
% Or = So(:,nx+1:end);
% Sf = Ox;

seed = 2233;

rng(seed)
xini = sign(randn(25,1));
% xini = 10*(randn(25,1));

