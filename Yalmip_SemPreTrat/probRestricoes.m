classdef probRestricoes <handle
properties 
    Robot
    rho
    rho_x
    P
    maxN
    Problem
    Sol
    SolTime
end
methods
    function obj = probRestricoes()
        
        obj.Robot = robot_model;
        obj.rho = [];
        obj.rho_x = [];
        obj.P = [];
        obj.maxN = [];
        obj.Problem = [];
        obj.Sol = [];
        obj.SolTime = 0;


    end
    function obj = create_yalmip_prob(obj)
        yalmip('clear');
        constraints = [];
        objective = 0;

        % 1. Inicializacao
        rbar_var = {sdpvar(1,1)};
        N_var = intvar(repmat(1,1,2),ones(1,2));
        U = sdpvar(repmat(3,1,obj.maxN),ones(1,obj.maxN));
        X = sdpvar(repmat(6,1,obj.maxN+1),ones(1,obj.maxN+1));
        
        rbar = [rbar_var{1}; rbar_var{1} + 2*pi*N_var{1}; rbar_var{1} + 2*pi*N_var{2}];

        xbar = obj.Robot.State_ref*rbar; ubar = obj.Robot.Input_ref*rbar;
        bf = obj.Robot.Term_bo - obj.Robot.Term_Or*rbar;

    N = obj.maxN;
    for k = 1:N

    objective = objective + (X{k} - xbar)'* obj.rho_x*(X{k} - xbar) + (U{k} - ubar)'*obj.rho*(U{k} - ubar);

    constraints = [constraints,-5 <= U{k} <= 5];
    constraints = [constraints,X{k+1} ==  obj.Robot.A*X{k}+obj.Robot.B*U{k}];
    end

    
    objective = objective + (X{N+1} - xbar)'*obj.P*(X{N+1} - xbar);
    constraints = [constraints, obj.Robot.Term_Sf*X{N+1} <= bf];

    ops = sdpsettings('Verbose',0,'solver','gurobi','gurobi.IntFeasTol',1e-8);
    obj.Problem = optimizer(constraints,objective,ops,X{1},[X U rbar_var N_var]);
    end

    function obj = solve(obj)
        X0 = cell(1,1);
        X0{1} = obj.Robot.x0;
        tic;
        Sol = obj.Problem(X0);
        obj.SolTime = toc;
        obj.Sol = Sol;
        X = [Sol{1:obj.maxN+1}];
        obj.Robot.setX(X);
        U = [Sol{(obj.maxN+1)+1: obj.maxN+1+obj.maxN}];
        obj.Robot.setU(U);
        rbar_var = [Sol{obj.maxN*2+2}];
        obj.Robot.setrbar_var(rbar_var);
        N_var = [Sol{obj.maxN*2+3:end}];
        obj.Robot.setN_var(N_var);
    end
end

end