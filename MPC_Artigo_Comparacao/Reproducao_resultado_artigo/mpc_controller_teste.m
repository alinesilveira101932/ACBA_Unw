function controller = mpc_controller_teste(A,B,N,Q,R,S,na,nx,nu,I_na,beta)

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
        constraints = [constraints, -0.3*ones(10,1) <= U{i}, U{i}  <= 0.3*ones(10,1)];
        constraints = [constraints, X{i+1} == kron(I_na,A)*X{i}+kron(I_na,B)*U{i}];
        objective = objective + (X{N+1})'*S*(X{N+1});
        
end

% Final
    
    constraints = [constraints, (X{N+1})'*S*(X{N+1})<= beta^2];

    
% Solucionar
ops = sdpsettings('verbose',0,'solver','gurobi');
outputs = [X U];
controller = optimizer(constraints,objective,ops,x0,outputs);

end

 % for i =1:100
 %     x0 = sign(randn(25,1));
 %     val(i) = max(abs(K*x0));
 % end