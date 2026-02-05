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
nr = 3; % numero de referencias

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

%% Restricoes 

% Su *u <= bu
% Ag1
umax1 = 1000;  umin1 = -1000;
% Ag2
umax2 = 1000;  umin2 = -1000;
% Ag3
umax3 = 10000;  umin3 = -1000;

% Sx*x <= bx
% Ag1
x1max_1 = 1000; x1min_1 = -1000;
x2max_1 = 1000; x2min_1 = -1000;
% Ag2
x1max_2 = 1000; x1min_2 = -1000;
x2max_2 = 1000; x2min_2 = -1000;
% Ag3
x1max_3 = 1000; x1min_3 = -1000;
x2max_3 = 1000; x2min_3 = -1000;

% Matrizes empilhando agentes
Su = [1; -1];
eu = 1e-3*ones(na*size(Su,1),1);

Su_barra = Su;
for i = 2:na
    Su_barra = blkdiag(Su_barra,Su);
end

Sx = [1 0; -1 0; 0 1; 0 -1]; 
ex = 1e-3*ones(na*size(Sx,1),1);

Sx_barra = Sx;
for i = 2:na
    Sx_barra = blkdiag(Sx_barra,Sx);
end

bx_barra = [x1max_1 ; -x1min_1; x2max_1; -x2min_1;
            x1max_2 ; -x1min_2; x2max_2; -x2min_2;
            x1max_3 ; -x1min_3; x2max_3; -x2min_3;
                                                    ];
bu_barra = [umax1; -umin1;
            umax2; -umin2;
            umax3; -umin3 ];