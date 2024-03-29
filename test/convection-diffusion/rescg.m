function [x, k, res, resvec] = rescg(A, b, restart, tol, x0, maxit)

%  RESCG   Restarted semi-conjugate gradient method
%
%   [x] = rescg(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b. The n-by-n matrix A must be positive definite but need not be symmetric.
%   The right-hand side column vector b must have length n.
%
%   [x] = rescg(A, b, restart) specifies the restart number. If restart is []
%   then rescg uses the default, n.  (This is usually too large when n is large.)
%
%   [x] = rescg(A, b, restart, tol) specifies the stopping tolerance.
%   If tol is [], rescg uses the default, 1e-6.
%   Note: Even if the user specifies x0, the stopping rule for Ax = b is relevant:
%         norm(r)/norm(b) <= tol, where r = b - A(x0+d) and d is called x below.
%
%   [x] = rescg(A, b, restart, tol, x0) specifies the initial guess.
%   If x0 is [], rescg uses the default, an all zero vector.
%
%   [x] = rescg(A, b, restart, tol, x0, maxit) specifies the iteration limit.
%   If maxit is [], rescg uses the default, 10000.
%
%   [x, k] = rescg(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = rescg(A, b, ...) returns the last relative residual norm
%   res = norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = rescg(A, b, ...) returns a vector of estimates of the
%   residual norms norm(b-Ax) at each iteration.
%------------------------------------------------------------------
% 2022-10-17
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


resvec = zeros(maxit, 1);


r   = b;
nr  = nb;
tol = tol*nb;
k   = 0;

%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    p  = r;
    np = norm(p);
    p  = p/np;
    q  = A*p;
    sp = p'*q; %%for reducing computation
    L  = sp;
    P  = p;
    Q  = q;
    for inner = 1:restart
        k = k+1;
        resvec(k) = nr;
        a  = (p'*r)/sp;
        x  = x + a*p;
        r  = r - a*q;
        nr = norm(r);
        if (nr < tol)
            break;
        end
        v   = A*r;
        lam = L\(P'*v);
        p   = r - P*lam;
        np  = norm(p);
        p   = p/np;
        q   = v - Q*lam;
        q   = q/np;
        sp  = p'*q; %for reducing computation
        L   = [L, sparse(inner,1); p'*Q, sp];
        P   = [P,p];
        Q   = [Q,q];
    end
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end