function [x, k, res, resvec] = swiwp(A, b, m, tol, x0, maxit)

%  SWIWP   Sliding window implementation without pre-allocated memory
% 
%   [x] = swiwp(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b. The n-by-n coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector b must have length n.
%
%   [x] = swiwp(A, b, m) specifies the number of the sliding window. If m is [] 
%   then swiwp uses the default, n.
%
%   [x] = swiwp(A, b, m, tol) specifies the tolerance of the method. If tol is []
%   then swiwp uses the default, 1e-6.
%
%   [x] = swiwp(A, b, m, tol, x0)  specifies the initial guess. If x0 is [] 
%   then swiwp uses the default, an all zero vector.
%
%   [x] = swiwp(A, b, m, tol, x0, maxit) specifies the maximum number of iterations.
%   If maxit is [] then swiwp uses the default, 10000.
%
%   [x, k] = swiwp(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = swiwp(A, b, ...) also returns the last relative
%   residual norm norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = swiwp(A, b, ...) also returns a vector of estimates of the 
%   residual norms at each iteration, including norm(b-Ax).
%------------------------------------------------------------------
% 2021-06-01
% Na Huang


n  = length(b);    
nb = norm(b);


%------------------------------------------------------------------
%     Retrieve input arguments.
%-----

if (nargin < 2)
   resvecor('MATLAB:NotEnoughInputs', 'Not enough input arguments.');
end
if (nargin < 3) || isempty(m)
   m = n;
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
x = zeros(n,1);


resvec = zeros(maxit, 1);
o      = sparse(m-1,1);


r    = b - A*x; 
nr   = norm(r);
p    = r;
np   = norm(p);
p    = p/np; 
q    = A*p;
sp   = p'*q; %%for reducing computation  
L    = sp;
P    = p; 
Q    = q; 
vm  = 2 : m;
tol  = tol*nb;

k    = 0;

%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    k = k+1;   
    resvec(k) = nr;
    % Compute normalized search direction. 
    a    = (p'*r)/sp;
    x    = x + a*p;
    r    = r - a*q;   
    nr   = norm(r);
    v    = A*r;
    lam  = L\(P'*v); 
    p  = r - P*lam; 
    np = norm(p); 
    p  = p/np;
    q  = v - Q*lam;  
    q  = q/np;
    sp = p'*q; %for reducing computation  
    if k < m
        L = [L, sparse(k,1); p'*Q, sp]; 
        P = [P,p]; 
        Q = [Q,q];
    else
        Pj = P(:,vm); 
        Qj = Q(:,vm);
        L  = [L(vm,vm), o; p'*Qj, sp]; 
        P  = [Pj,p]; 
        Q  = [Qj,q];
    end
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end