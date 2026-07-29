% Parametros dados - artigo - simulacao
clc; clear;
S2 = [2.551,-0.447,0.119,-0.813,-1.069;
     -0.447,4.028,0.227,1.356,-2.664;
     0.119,0.227,1.799,0.740,-2.431;
     -0.813,1.356,0.740,3.884,-3.689;
     -1.069,-2.664,-2.431,-3.689,10.081];

Q2 = [1,0,0,0,-1;
     0,1,0,0,-1;
     0,0,1,0,-1;
     0,0,0,1,-1;
    -1,-1,-1,-1,4];

alpha = 10; c = 0.152 ; W = 0.5*eye(5);

N = 9;  beta = 28.8;
T = 0.1;

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

% -0.3 ⩽u1 ⩽ 0.3, and -0.3 ⩽ u2 ⩽ 0.3. 

 L = [2,-1,0,-1,0;
     -1,2,-1,0,0
     0,-1,2,0,-1;
     -1,0,0,2,-1;
     0,0,-1,-1,2];

  % Dimencoes da planta
nx = length(A);
nu = size(B,2);

na = 5; % numero de agentes
I_na = eye(na);

 %% Calculos
 S1 = W*L;
 S = kron(S1,S2);
 R2 = alpha*B'*S2*B;
 K2 = -inv(B'*S2*B+R2)*B'*S2*A;
 K = c*kron(L,K2);
 
 Z = null(S);

 for i = 1:5
    nullK(:,i)= K*Z(:,i);
 end

% Theorem 8
R1 = W*((1+alpha)*I_na-c*L)/(c*alpha);
H = A'*S2*B*inv(B'*S2*B)*B'*S2*A;
R = kron(R1,R2);
pond = (c*S1*L)/(1+alpha);
Q = kron(S1,Q2)+kron(pond,H);
 
 % a = max(nullK(:))
 % b= min(nullK(:))

umax = 0.3;  umin = -0.3;
Su = [1 0; -1 0;0 1;0 -1];
b_u = [umax; -umin;umax;-umin];
Au = kron(I_na,Su);
bu = kron(ones(na,1),b_u);

%% Determinar O_inf
% Formula do maior elipsoide inscrito em um poliedro
% Boyd - Convex Optimization p. 414
V = orth(S');
% x = V*z -> z'V'S V z <= beta
% K V z pertence a U
% 
Kv = K*V;
Sv = V'*S*V;

% for i = 1:2*na
% beta(i,:) = bu(i,:)^2/(Au(i,:)*Kv*inv(Sv)*Kv'*Au(i,:)');
% end
% betai = beta(isfinite(beta));
% beta = min(betai(:));
% beta = 2;
% 
seed = 2233;

rng(seed)
x(:,1) = sign(randn(25,1));