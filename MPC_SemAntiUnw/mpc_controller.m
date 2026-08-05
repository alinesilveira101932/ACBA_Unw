function controller = mpc_controller(A,B,N,Q,R,P,W,bo,Oxbar,Oubar,Sf,na,nx,nu,Su,bu)

yalmip('clear');

x0 = sdpvar(na*nx,1);
xbar_old = sdpvar(nx,1);

% Variaveis
% Parametrizado em termos de q para igualdade de equilibrio no conjunto
% terminal
% qbar = sdpvar(nt,1);
% xbar = Tx*qbar;
% ubar = Tu*qbar;
% 
% % Equilibrios
% bf = bo - Oq*qbar;

% Variaveis
xbar = sdpvar(nx,1);
ubar = sdpvar(nu,1);

% Equilibrios
bf = bo - Oxbar*xbar-Oubar*ubar;

% Variaveis
for i = 1: na
  U{i} = sdpvar(nu,N);
  X{i} = sdpvar(nx,N+1);
end
constraints = [];
objective = 0;

% Inicial
for i = 1:na
constraints = [constraints, X{i}(:,1) == x0((i-1)*nx+1:i*nx)];
end

% Ao longo do horizonte
for i = 1:N
    for j = 1:na
        objective = objective + (X{j}(:,i)-xbar)'*Q*(X{j}(:,i)-xbar) + (U{j}(i)-ubar)'*R*(U{j}(i)-ubar);
        constraints = [constraints, -Su*U{j}(i) <= bu];
        constraints = [constraints, X{j}(:,i+1) == A*X{j}(:,i)+B*U{j}(i)];
        constraints = [constraints, (A-eye(nx))*xbar + B*ubar == 0];
    end
end

% Final
for j = 1:na
    objective = objective + (X{j}(:,N+1)-xbar)'*P*(X{j}(:,N+1)-xbar);
    constraints = [constraints, Sf*X{j}(:,N+1) <= bf];
end
objective = objective + (xbar-xbar_old)'*W*(xbar-xbar_old);

% Solucionar
ops = sdpsettings('verbose',0,'solver','gurobi');
outputs = [X U {xbar} {ubar}];
controller = optimizer(constraints,objective,ops,[x0; xbar_old],outputs);

end