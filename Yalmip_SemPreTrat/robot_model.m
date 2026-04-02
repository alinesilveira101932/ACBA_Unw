classdef robot_model < handle
% Robot model used in the optmization
properties
    A   % State matrix
    B   % Input matrix
    State_ref   % states variables' reference
    Input_ref   % input variables' reference
    Term_bo
    Term_Or 
    Term_Sf
    x0          % current state
    X           % predicted trajectory
    U           % predicted input sequence
    rbar_var
    N_var
end
methods
    function obj = robot_model()
        obj.A = [];
        obj.B = [];
        obj.Input_ref = [];
        obj.State_ref = [];
        obj.Term_bo = [];
        obj.Term_Or = [];
        obj.Term_Sf = [];
        obj.x0 = [];
        obj.X = [];
        obj.U = [];
        obj.rbar_var = [];
        obj.N_var = [];
    end
    function obj = evolve_state(obj,u)
        x = obj.A*obj.x0 + obj.B*u;
        obj.x0 = x;
    end
    function obj = setX(obj,X)
        obj.X = X;
    end
    function obj = setU(obj,U)
        obj.U = U;
    end
    function obj = setrbar_var(obj,rbar_var)
        obj.rbar_var = rbar_var;
    end
    function obj = setN_var(obj,N_var)
        obj.N_var = N_var;
    end
end

end