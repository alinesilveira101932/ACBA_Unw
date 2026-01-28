%% Parametros - Artigo MPC + Unwinding
% Auxiliar
% Autor: Aline Isabel
% Data: 23/01/2026

%--------------------------------------------------------------------
% Planta
Ac = [0 1; 0 0];
Bc = [0;1];
Cc = [1 0];

T = 0.1; % Periodo de amostragem

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
nq = size(Cc,1);

na = 3; % numero de agentes

N = 25; % Horizonte de predicao

% Matrizes de Peso
Q = 10*eye(nx);
R = 1;

% Discretizacao
[Phi,Gamma] = c2dm(Ac,Bc,[],[],T, 'zoh');

% Empilhando as matrizes para todos os agentes
A = kron(eye(na),Phi);
B = kron(eye(na),Gamma);
C = kron(eye(na),Cc);

Ac_barra = kron(eye(na),Ac);
Bc_barra = kron(eye(na),Bc);

R_barra = R;
for i = 2:na
    R_barra = blkdiag(R_barra,R);
end

Q_barra = Q;
for i = 2:na
    Q_barra = blkdiag(Q_barra,Q);
end

% LQR
[K,P_barra] = dlqr(A,B,Q_barra,R_barra);
Af = A-B*K;

% Valores de equilibrio
aux1 = inv([A-eye(nx*na) B;C zeros(na)])*[zeros(nx*na,na);eye(na)];
Nx = aux1(1:nx*na,:); Nu = aux1(nx*na+1:end,:);
Nx = Nx*[1;1;1]; % Mapeia pra uma unica ref
Nu = Nu*[1;1;1]; % Mapeia pra uma unica ref
%% Restricoes 
% Su *u <= bu 
umax = 10000;
umin = -1000;

Su = [1; -1];
bu = [umax; -umin];
eu = 1e-3*ones(na*size(Su,1),1);

% Sx*x <= bx
x1max = 1000;
x1min = -1000;
x2max = 1000;
x2min = -1000;

Sx = [1 0; -1 0; 0 1; 0 -1]; 
bx = [x1max ; -x1min; x2max; -x2min];
ex = 1e-3*ones(na*size(Sx,1),1);

% Matrizes empilhando agentes
Su_barra = Su;
for i = 2:na
    Su_barra = blkdiag(Su_barra,Su);
end

Sx_barra = Sx;
for i = 2:na
    Sx_barra = blkdiag(Sx_barra,Sx);
end

bx_barra = [repmat(bx,na,1)];
bu_barra = [repmat(bu,na,1)];