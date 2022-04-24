function [x, k, res, resvec] = swiwp(A, b, mk, tol, x0, maxit)
%% 2021-06-01
%% Na Huang

%  SWI   Sliding window implementation without pre-allocated memory
% 
%   [Z] =  SWIWP(A, B) attempts to find a solution Z to the system of linear equations
%   AZ=B.  The N-by-N coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector B must have length N.
%
%   [Z] =  SWIWP(A, B, MK) specifies the number of the sliding window. If MK is []
%   then  SWIWP uses the default, N.
%
%   [Z] =  SWIWP(A, B, MK, TOL) specifies the tolerance of the method. If TOL is []
%   then  SWIWP uses the default, 1e-6.
%
%   [Z] =  SWIWP(A, B, MK, TOL, Z)  specifies the initial guess.  If Z is [] then 
%    SWIWP uses the default, an all zero vector.
%
%   [Z] =  SWIWP(A, B, MK, TOL, Z, MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then  SWIWP uses the default, 10000.
%
%   [Z, K]= SWIWP(A, B, ...) returns the iteration number at which Z
%   was computed: 1 <= ITER <= MAXIT.
%
%   [Z, K, RES]= SWIWP(A, B, ...) also returns the last relative
%   residual norm NORM(B-AZ)/NORM(B).
%
%   [Z, K, RESVEC]= SWIWP(A, B, ...) also returns a vector of estimates of the 
%   residual norms at each iteration, including NORM(B-AZ).
%------------------------------------------------------------------


n  = length(b);    
nb = norm(b);


%------------------------------------------------------------------
%     Retrieve input arguments.
%-----

if (nargin < 2)
   resvecor('MATLAB:NotEnoughInputs', 'Not enough input arguments.');
end
if (nargin < 3) || isempty(mk)
   mk = n;
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

%% sovle Ad = A(x-x0) = b-Ax0 = r0 first
if gotx0 ~= 0
    b = b - A*x0;
end
x = zeros(n,1);


resvec = zeros(maxit, 1);
o      = sparse(mk-1,1);


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
vmk  = 2 : mk;
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
    sp = p'*q; %%for reducing computation  
    if k < mk
        L = [L, sparse(k,1); p'*Q, sp]; 
        P = [P,p]; 
        Q = [Q,q];
    else
        Pj = P(:,vmk); 
        Qj = Q(:,vmk);
        L  = [L(vmk,vmk), o; p'*Qj, sp]; 
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