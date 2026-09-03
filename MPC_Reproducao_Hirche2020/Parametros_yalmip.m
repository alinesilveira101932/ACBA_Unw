%% Parametros - Usando YALMIP - Artigo MPC + Unwinding - Agentes separados
% Auxiliar
% Autor: Aline Isabel
% Data: 31/03/2026

%--------------------------------------------------------------------
% Planta
Ac = [0 1; 0 0];
Bc = [0;1];
C = [1 0];

% % b = 0.5;
% % Ac = [0 1; 0 -b];
% Bc = [0;1];
% C = [1 0];

T = 0.1; % Periodo de amostragem

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
ny = size(C,1);

na = 3; % numero de agentes

N = 40; % Horizonte de predicao

% Matrizes de Peso
Q = 0.4*[5 0; 0 1];
R = 0.02;

% Discretizacao
[A,B] = c2dm(Ac,Bc,[],[],T, 'zoh');

% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;

%% Restricoes 

% Su *u <= bu
umax = 5;  umin = -5;

% Sx*x <= bx
x1max = 10; x1min = -10;
x2max = 10; x2min = -10;

Su = [1; -1];
eu = 1e-3*ones(1*size(Su,1),1);

Sx = [1 0; -1 0; 0 1; 0 -1]; 
ex = 1e-3*ones(1*size(Sx,1),1);

bx = [x1max ; -x1min; x2max; -x2min];
bu = [umax; -umin];

%% Determinar O_psi_inf - com respeito a q -> xbar e ubar

H = [A-eye(nx) B]; T = null(H);
Tx = T(1:nx,:); Tu = T(nx+1:end,:);
nt = size(T,2);

Gamma = [eye(nx) zeros(nx,1); -K K*Tx+Tu; zeros(nx,nx) Tx; zeros(nu,nx) Tu];

Spsi = blkdiag(Sx,Su,Sx,Su);
bpsi = [bx; bu; bx-ex; bu-eu];

Apsi_f = [Af B*K*Tx+B*Tu; zeros(nu,nx) eye(nu)];
max_iter = 180;
tol = 1e-3;

[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter,tol);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

% figure (1)
% plot(O_psi_inf)

% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Oq = So(:,nx+1:end);
Sf = Ox;

%% Conjunto Y
for i = 1:na
% Projecao de Oinf em q_i:
S_q = O_psi_inf.projection(nx+1 : nx+nt);
S_q.minHRep();

% Mapeia para o espaço de saida: ybar = C*xs = C*Tx*q
M = C*Tx;
Omega_i_s = S_q.affineMap(M);   
end

% Intersecção global
Ys = Omega_i_s;
for i = 2:na
    Ys = Ys & Omega_i_s;
end
Ys.minHRep();   
Ay = Ys.A;              
by = Ys.b;             
