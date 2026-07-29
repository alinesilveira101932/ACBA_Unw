function controller = mpc_controller(A,B,N,Q,R,S,na,nx,nu,Au,bu,I_na,beta)

yalmip('clear');

x0 = sdpvar(na*nx,1);

% Variaveis
  U = sdpvar(repmat(na*nu,1,N),ones(1,N));
  X = sdpvar(repmat(na*nx,1,N+1),ones(1,N+1));

constraints = [];
objective = 0;

% Inicial
constraints = [constraints, X{1} == x0];

% Ao longo do horizonte
for i = 1:N
        objective = objective + X{i}'*Q*X{i} + U{i}'*R*U{i};
        constraints = [constraints, Au*U{i} <= bu];
        constraints = [constraints, X{i+1} == kron(I_na,A)*X{i}+kron(I_na,B)*U{i}];
        
end

% Final
    
    constraints = [constraints, (X{N+1})'*S*(X{N+1})<= beta^2];
    objective = objective + (X{N+1})'*S*(X{N+1});
% Solucionar
ops = sdpsettings('verbose',0,'solver','gurobi');
outputs = [X U];
controller = optimizer(constraints,objective,ops,x0,outputs);

end