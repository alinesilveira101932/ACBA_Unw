%% Matrizes de Parametros - Artigo MPC + Unwinding
% Auxiliar
% Autor: Aline Isabel
% Data: 23/01/2026

%% Ap
Ap = A;
for i = 2:N
    Ap = [Ap;A^i];
end

%% Bp

for i = 1:N
    for j = 1:N
        if i >= j
            Bp{i,j} = A^(i-j)*B;
        else
            Bp{i,j} = zeros(na*nx,1*na);
        end
    end
end
Bp = cell2mat(Bp);


%% Matrizes de peso
Rp = R_barra;
for i = 2:N
    Rp = blkdiag(Rp,R_barra);
end

Qp = Q_barra;
for i = 2:N-1
    Qp = blkdiag(Qp,Q_barra);
end
Qp = blkdiag(Qp,P_barra); % P é Restricao final para k = N

%% Matrizes das restricoes de U
% Sup * u < bup

Sup = Su_barra;
for i = 2:N  
    Sup = blkdiag(Sup,Su_barra);
end

bup = repmat(bu_barra,N,1);

%% Matrizes das restricoes de X
% Sxp * x < bxp
Sxp = Sx_barra;
for i = 2:N-1
    Sxp = blkdiag(Sxp,Sx_barra);
end

Sxp = blkdiag(Sxp,Sf); % Sf*x < bf é Restricao final para k = N

%% Auxiliar
Inx = repmat(eye(nx*na),N,1);
Inu = repmat(eye(nu*na),N,1);
Nu_barra = repmat(Nu,N,1);
Nx_barra = repmat(Nx,N,1);
