function [x, k, res, resvec] = refom(A, b, restart, tol, x0, maxit)

% refom Restarted Full Orthogonalization Method
% Correspond to Algorithms 6.4 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"
% solve the triangular system used backward substitution

%   [x] = refom(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b. The n-by-n coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector b must have length n.
%
%   [x] = refom(A, b, tol) specifies the tolerance of the method. If tol is []
%   then refom uses the default, 1e-6.
%
%   [x] = refom(A, b, tol, x0)  specifies the initial guess. If x0 is []
%   then refom uses the default, an all zero vector.
%
%   [x] = refom(A, b, tol, x0, maxit) specifies the maximum number of iterations.
%   If maxit is [] then refom uses the default, 10000.
%
%   [x, k] = refom(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = refom(A, b, ...) also returns the last relative
%   residual norm norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = refom(A, b, ...) also returns a vector of estimates of the
%   residual norms at each iteration, including norm(b-Ax).
%------------------------------------------------------------------
% 14 Sep 2022
% Na Huang

n  = length(b);


%------------------------------------------------------------------
%     Retrieve input arguments.
%-----

if (nargin < 2)
    resvecor('MATLAB:NotEnoughInputs', 'Not enough input arguments.');
end
if (nargin < 3) || isempty(restart)
    restart = n;
end
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
if (nargin < 5) || isempty(x0)
    gotx0 = 0;
end
if (nargin < 6) || isempty(maxit)
    maxit = 10000;
end

%------------------------------------------------------------------
%     Initialization.
%-----

% sovle Ad = A(x-x0) = b-Ax0 = r0 first
if gotx0 ~= 0
    b = b - A*x0;
end
nb = norm(b);
x  = zeros(n,1);

%------------------------------------------------------------------
%     Initialization.
%-----



resvec = zeros(maxit, 1);
e      = zeros(restart , 1);
r      = b;
nr     = nb;
tol    = tol*nb;
k      = 0;


%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    beta   = nr;
    V(:,1) = r/beta;
    e(1)   = beta;
    for inner = 1:restart
        k = k+1;
        resvec(k) = nr;
        %% Orthogonalization Process
        w = A*V(:,inner);
        for j = 1:inner
            Vj = V(:,j);
            H(j,inner) = w'*Vj;
            w = w - H(j,inner)*Vj;
        end
        nw = norm(w);
        H(inner+1,inner) = nw;
        V(:,inner+1) = w/H(inner+1,inner);
        %% Use Givens transformations to transform H(1:inner, 1:inner) into an upper triangular
        Hkk = H(1:inner,1:inner);
        ge  = e(1:inner);
        for i = 1:(inner-1)
            veci = i:(i+1);
            G  = givens(Hkk(i,i),Hkk(i+1,i));
            Hkk(veci,:) = G*Hkk(veci,:);
            ge(veci) = G*ge(veci);
        end
        ykinner = ge(inner) / Hkk(inner, inner);
        nr = nw*abs(ykinner);   % Proposition 6.7 & formula (6.16)
        if (nr < tol)
            break;
        end
    end
    %     %% yk = Hkk\ge;  Use bacinnerward substitution to solve the upper triangular system
    %     yk(inner) = ykinner;
    %     for i = (inner-1) : -1 : 1
    %         vi = i+1 : inner;
    %         syi = Hkk(i, vi) * yk(vi);
    %         yk(i) = (ge(i)-syi)/Hkk(i,i);
    %     end
    yk = sparse(Hkk)\ge;
    x  = x + V(:,1:inner)*yk;
    r  = b - A*x;
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end