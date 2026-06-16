%% Parametros - DOI: 10.1561/2600000019.
% Auxiliar
% Autor: Aline Isabel
% Data: 05/2026
clc;clear
%--------------------------------------------------------------------
%% Planta
% Integrador duplo - sistema instavel
% Ac = [0 1; 0 0];
% Bc = [0;1];
% C = [1 0];

% Integrador duplo com atrito viscoso - sistema semi-estavel
b = 0.5;
Ac = [0 1; 0 -b];
Bc = [0;1];
C = [1 0];

T = 0.1; % Periodo de amostragem

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
nq = size(C,1);

na = 3; % numero de agentes
nr = 3; % numero de referencias

N = 25; % Horizonte de predicao

% Discretizacao
[A,B] = c2dm(Ac,Bc,[],[],T, 'zoh');

% Sistema empilhado
I_na = eye(na);

sistema = 0; % 0 para sistemas semi-estaveis e 1 para sistemas instaveis - definicao artigo Li2018

%% Matrizes do grafo - EqsRef:  DOI: 10.1561/2600000019.
%        1
%        o
%       / \
%      o - o
%      2   3
% todos se comunicam

% The adjacency matrix of a digraph is defined by 
% Adj =[a_{ij}] ∈ R^{n×n} where a_{ij} = 1 if(j,i) ∈ E, and a_{ij} = 0 otherwise.
Adj = [0 1 1 
       1 0 1
       1 1 0]; 

% The in-degree matrix of a digraph is defined as
% Din diag([din_1 ,...,din_n]), where din_i = \sum_(j=1)^(n) a_{ij}
Din = [2 0 0
       0 2 0
       0 0 2];

% The Laplacian matrix of the digraph is then defined by L = Din−Adj,
L =  [2 -1 -1
      -1 2 -1
      -1 -1 2]; 

%% Sistemas instaveis
if sistema == 1 

% Matrizes de Peso - Conjunto PB do artigo (instavel) - dois polos em 1
Q2 = eye(nx); %(A,Q2) observavel
alpha = 0.5; % Constante alpha > 0 exigida pelo Lemma 7

W = 2*eye(na);
S1 = W*L;

% Para o integrador duplo,lambda = 1 -> delta_c = 1 - 1/(max|lambda_u|^2) = 1 - 1/1 = 0.
delta = 0.6; % delta > delta_c = 0

% Resolve a ARE do lemma 7
% R_pond = ((1+alpha)/delta) -1; % deve ser somado ao termo do meio de forma que apareca a ponderacao desejada no lemma
% 
% [S2,~,~]=idare(A,B,Q2,R_pond,[],[]);
rho = delta/(1+alpha);
E = eye(2);
S2 = eye(nx);

for k=1:100

    Rk = (1/rho - 1)*(B'*S2*B);
    Snew = idare(A,B,Q2,Rk,[],E);

    if norm(Snew-S2,'fro') < 1e-10
        break
    end

    S2 = Snew;
end
Residual =A'*S2*A - S2 + Q2 - (delta/(1+alpha))*A'*S2*B*inv(B'*S2*B)*B'*S2*A;
norm(Residual);

% delta/min(eig(L) < c < 1/max(eig(L)  -> eig(L) = [0 3 3]
% 0,2 < c < 0,33
c = 0.3;

R1 = W*((1+alpha)*eye(na)-c*L)/(c*alpha);
R2 = alpha*B'*S2*B;

% Theorem 8
H = A'*S2*B*inv(B'*S2*B)*B'*S2*A;
R = kron(R1,R2);
pond = (c*S1*L-delta*S1)/(1+alpha);
Q = kron(S1,Q2)+kron(pond,H);

K2 = -inv(B'*S2*B+R2)*B'*S2*A;

K = c*kron(L,K2);
S = kron(S1,S2);

%% Sistemas semi-estaveis
elseif sistema ==0
C2 = [0 0.1]; % (A,C2) semi observavel -> rank = n-1
Q2 = C2'*C2;
% Q2 = eye(2);
% Lemma 6: A'SA− S = C'C using X = dlyap(A,Q) solves the discrete-time Lyapunov equation AXA' − X + Q = 0,
s11 = 1; % existem varias solucoes
syms s12 s22
S2 = [s11 s12; s12 s22];
eq = simplify(A'*S2*A-S2+Q2);
sol = solve(eq(1,1)==0,eq(1,2)==0,eq(2,2)==0);
s12 = double(sol.s12); s22 = double(sol.s22);

S2 = [s11 s12; s12 s22]
eig(S2)

W = eye(na);
S1 = W*L;
% < c < 1/max(eig(L)  -> eig(L) = [0 3 3] -> c < 0,33
c = 0.3;
alpha = 1; % Constante alpha > 0 
R1 = W*((1+alpha)*I_na-c*L)/(c*alpha);
R2 = alpha*B'*S2*B;

% Teorema 8
S = kron(S1,S2);
K2 = -inv(B'*S2*B+R2)*B'*S2*A;
K = c*kron(L,K2);
H = A'*S2*B*inv(B'*S2*B)*B'*S2*A;
R = kron(R1,R2);
pond = (c*S1*L)/(1+alpha);
Q = kron(S1,Q2)+kron(pond,H);

Residual = A'*S2*A - S2 + Q2;
norm(Residual)

end

%% Restricoes 
% Su_i *u <= bu_i
% [1; -1] u <= [umax; -umin];

umax = 5;  umin = -5;
Su = [1; -1];
b_u = [umax; -umin];
Au = kron(I_na,Su);
bu = kron(ones(na,1),b_u);

%% Determinar O_inf
% Formula do maior elipsoide inscrito em um poliedro
% Boyd - Convex Optimization p. 414
V = orth(S');
% x = V*z -> z'V'S V z <= beta
% K V z pertence a U

Kv = K*V;
Sv = V'*S*V;

for i = 1:2*na
beta(i,:) = bu(i,:)^2/(Au(i,:)*Kv*inv(Sv)*Kv'*Au(i,:)');
end
betai = beta(isfinite(beta));
beta = min(betai(:));
