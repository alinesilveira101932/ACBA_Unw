function controller = mpc_controller(A,B,N,Q,R,P,W,bo,Oxbar,Oubar,Sf,na,nx,nu,Su,bu)

yalmip('clear');

x0 = sdpvar(na*nx,1);

% Variaveis
xbar  = sdpvar(nx,1); % comum a todos

for i = 1: na
  U{i} = sdpvar(nu,N);
  X{i} = sdpvar(nx,N+1);
  XBAR_old{i} = sdpvar(nx,1); % pra cada agente unw
  UBAR{i} = sdpvar(nu,1); % pra cada agente unw 
end

XBAR{1} = xbar; %definicao
for i = 1: na-1
   Nbar{i} = [intvar(1,1); 0]; % Unw apenas em x1 % pra cada agente unw
   XBAR{i+1} = xbar + 2*pi*Nbar{i}; % pra cada agente unw
end

constraints = [];
objective = 0;

% Inicial
for i = 1:na
constraints = [constraints, X{i}(:,1) == x0(nx*i-1:nx*i)];
end

% Equilibrios
for j = 1:na
    bf{j}   = bo - Oxbar*XBAR{j}-Oubar*UBAR{j}; % pra cada agente unw
end

% Ao longo do horizonte
for i = 1:N
    for j = 1:na
        objective = objective + (X{j}(:,i)-XBAR{j})'*Q*(X{j}(:,i)-XBAR{j}) + (U{j}(:,i)-UBAR{j})'*R*(U{j}(:,i)-UBAR{j});
        constraints = [constraints, Su*U{j}(:,i) <=bu];
        constraints = [constraints, X{j}(:,i+1) == A*X{j}(:,i)+B*U{j}(:,i)];
        constraints = [constraints, (A-eye(nx))*XBAR{j} + B*UBAR{j} == 0];
    end
end

% Final
for j = 1:na
    objective = objective + (X{j}(:,N+1)-XBAR{j})'*P*(X{j}(:,N+1)-XBAR{j});
    constraints = [constraints, Sf*X{j}(:,N+1) <= bf{j}];
    objective = objective + (XBAR{j}- XBAR_old{j})'*W*(XBAR{j}- XBAR_old{j});
end

% Solucionar
ops = sdpsettings('verbose',0,'solver','gurobi');
outputs = [X U XBAR UBAR Nbar];
ini = [x0];
for i = 1:na
    ini = [ini; XBAR_old{i}];
end
controller = optimizer(constraints,objective,ops,ini,outputs);

end