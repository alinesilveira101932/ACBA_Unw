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

N = 45; 

% A = 0.99*A;
% Dimencoes da planta
nx = length(A);
nu = size(B,2);
ny = size(C,1);

na = 5; % numero de agentes

% Matrizes de Peso
Q = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
W = 0.01;
R = [2.5*12.5144    5.8500
    5.8500   17.6198];
% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;

%% Restricoes 
umax = 0.30;  umin = -0.30;
Su = [1 0; -1 0;0 1;0 -1];
bu = [umax; -umin;umax;-umin];
eu = 1e-3*ones(1*size(bu,1),1);

Sxbar = [A-eye(nx);-A+eye(nx)];
bxubar = [zeros(nx,1); zeros(nx,1)];
Subar = [B;-B];

%% Determinar O_psi_inf - com respeito a q - xbar e ubar

H = [A-eye(nx) B]; T = null(H);
Tx = T(1:nx,:); Tu = T(nx+1:end,:);
nt = size(T,2);

Gamma = [-K K*Tx+Tu; zeros(nu,nx) Tu];

Spsi = blkdiag(Su,Su);
bpsi = [bu; bu-eu];

Apsi_f = [Af B*K*Tx+B*Tu; zeros(nu,nx) eye(nu)];
max_iter = 50;
tol = 1e-6;

[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter,tol);

O_psi_inf_1 = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox
% Operação de projeção empregando o MPT Toolbox

% D1 = projection(O_psi_inf_1,2:3);
% 
% % Apresentação gráfica de D e Oinf usando o plot do MPT Toolbox
% figure (89) ; plot(D1,'Color','g')


% figure (1)
% plot(O_psi_inf)

% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Oq = So(:,nx+1:end);
Sf = Ox;

% % Projecao Cond Inicial
% Sw = [(A-eye(nx)) B; -(A-eye(nx)) -B;zeros(nu*2,nx) Su];
% bw = [zeros(nx,1);zeros(nx,1);bu];
% 
% % Operação de projeção empregando o MPT Toolbox
% Pw = Polyhedron('H',[Sw bw]);
% Du = projection(Pw,nx+1:nx+nu);
% Dx = projection(Pw,1:nx);
% % Apresentação gráfica de D e Oinf usando o plot do MPT Toolbox
% figure (98); plot(Du,'Color','g')
% 
% 
% ubar01 = [-0.3; -0.2609]; ubar02 = [0.3; 0.2609];ubar03 = [0 ;0];ubar04 = [-0.3; -0.2609];ubar05 = [0.0575; 0.05];
% xbar01 = sdpvar(nx,1); xbar02 = sdpvar(nx,1);xbar03 = sdpvar(nx,1);xbar04 = sdpvar(nx,1);xbar05 = sdpvar(nx,1);
% F1 = [(A-eye(nx))*xbar01 == -B*ubar01];F2 = [(A-eye(nx))*xbar02 == -B*ubar02];F3 = [(A-eye(nx))*xbar03 == -B*ubar03];F4 = [(A-eye(nx))*xbar04 == -B*ubar04];F5 = [(A-eye(nx))*xbar05 == -B*ubar05];
% J1 = xbar01'*xbar01;J2 = xbar02'*xbar02;J3 = xbar03'*xbar03;J4 = xbar04'*xbar04;J5 = xbar05'*xbar05;
% 
% optimize(F1,J1);optimize(F2,J2);optimize(F3,J3);optimize(F4,J4);optimize(F5,J5);
% 
% xbar01 = value(xbar01);xbar02 = value(xbar02);xbar03 = value(xbar03);xbar04 = value(xbar04);xbar05 = value(xbar05);


%% Determinar O_psi_inf - com respeito a xbar -> sendo ubar = 0;

Gamma = [-K K eye(nu); zeros(nu,nx) zeros(nu,nx) eye(nu);zeros(nx,nx) eye(nx) zeros(nx,nu)];
Spsi_1 = blkdiag(Su,Su);
bpsi_1 = [bu;bu-eu];
Spsi = [Spsi_1 zeros(2*nu*nu,nx); [zeros(nu*nx,nu) Subar Sxbar]];
bpsi = [bu;bu-eu;bxubar];
Apsi_f = [Af B*K B;  zeros(nu,nx) zeros(nu,nx) eye(nu);zeros(nx,nx) eye(nx) zeros(nx,nu)];
max_iter = 100;
tol = 1e-8;
[So_2,bo_2,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter,tol);

O_psi_inf_2 = Polyhedron('H',[So_2 bo_2]); % Declaração de Oinf como objeto do MPT Toolbox

% D2 = projection(O_psi_inf_2,2:3);
% figure (90); plot(D2,'Color','g')
% figure (1)
% plot(O_psi_inf)

% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Oubar = So(:,nx+1:nx+nu);
Oxbar = So(:,nx+nu+1:end);
Sf = Ox;

%% Condicao inicial
seed = 2233;

rng(seed)

% xini = 10*(randn(25,1));
% xini =[xbar01;xbar02;xbar03;xbar04;xbar05]
xini = sign(randn(25,1));