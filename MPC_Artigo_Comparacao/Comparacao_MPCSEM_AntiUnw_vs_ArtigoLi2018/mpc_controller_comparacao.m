%% Otimizacao - codigo auxiliar
% MPC - Comparacao sem anti Unwinding  vs artigo li2018
% Modelo artigo Li2018 - 5 estados e 2 entradas
% Autor: Aline Isabel
% Data: 24/07/2026
%--------------------------------------------------------------------
function controller = mpc_controller_comparacao(A,B,N,Q,R,P,W,Nx,Nu,bo,Or,Sf,na,nx,nu,Su,bu)

yalmip('clear');

x0 = sdpvar(na*nx,1);
rbar_old = sdpvar(2,1);

% Variaveis
rbar  = sdpvar(2,1);

for i = 1: na
  U{i} = sdpvar(nu,N);
  X{i} = sdpvar(nx,N+1);
end
constraints = [];
objective = 0;

% Inicial
constraints = [constraints, X{1}(:,1) == x0(1:nx)];
constraints = [constraints, X{2}(:,1) == x0(nx+1:2*nx)];
constraints = [constraints, X{3}(:,1) == x0(2*nx+1:3*nx)];
constraints = [constraints, X{4}(:,1) == x0(3*nx+1:4*nx)];
constraints = [constraints, X{5}(:,1) == x0(4*nx+1:5*nx)];

% Equilibrios
for j = 1:na
    xbar{j} = Nx*rbar;
    % ubar{j} = Nu*rbar;
    ubar{j} = [0 0;0 0]*rbar;
    bf{j}   = bo - Or*rbar;
end

% Ao longo do horizonte
for i = 1:N
    for j = 1:na
        objective = objective + (X{j}(:,i)-xbar{j})'*Q*(X{j}(:,i)-xbar{j}) + (U{j}(:,i)-ubar{j})'*R*(U{j}(:,i)-ubar{j})+((rbar-rbar_old)'*W*(rbar-rbar_old));
        constraints = [constraints, Su*U{j}(:,i) <= bu];
        constraints = [constraints, X{j}(:,i+1) == A*X{j}(:,i)+B*U{j}(:,i)];
    end
end

% Final
for j = 1:na
    objective = objective + (X{j}(:,N+1)-xbar{j})'*P*(X{j}(:,N+1)-xbar{j});
    constraints = [constraints, Sf*X{j}(:,N+1) <= bf{j}];
end


% Solucionar
ops = sdpsettings('verbose',0,'solver','gurobi');
outputs = [X U {rbar}];
controller = optimizer(constraints,objective,ops,[x0; rbar_old],outputs);

end