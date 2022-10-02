function [x, k, res, resvec] = swi(A, b, m, tol, x0, maxit)

%  SWI   Sliding window implementation with pre-allocated memory
% 
%   [x] = swi(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b.  The n-by-n coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector b must have length n.
%
%   [x] = swi(A, b, m) specifies the number of the sliding window. If m is [] 
%   then swi uses the default, n.
%
%   [x] = swi(A, b, m, tol) specifies the tolerance of the method. If tol is []
%   then swi uses the default, 1e-6.
%
%   [x] = swi(A, b, m, tol, x0)  specifies the initial guess. If x0 is [] 
%   then swi uses the default, an all zero vector.
%
%   [x] = swi(A, b, m, tol, x0, maxit) specifies the maximum number of iterations.
%   If maxit is [] then swi uses the default, 10000.
%
%   [x, k] = swi(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = swi(A, b, ...) also returns the last relative
%   residual norm norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = swi(A, b, ...) also returns a vector of estimates of the 
%   residual norms at each iteration, including norm(b-Ax).
%------------------------------------------------------------------
% 2021-06-01
% Na Huang


n  = length(b);    



%------------------------------------------------------------------
%     Retrieve input arguments.
%-----

if (nargin < 2)
   error('MATLAB:NotEnoughInputs', 'Not enough input arguments.');
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
nb = norm(b);
x  = zeros(n,1);


L      = zeros(m,m);
P      = zeros(n,m);
Q      = zeros(n,m); 
resvec = zeros(maxit,1);



k   = 0;
r   = b;  
nr  = nb;
tol = tol*nb;
 

%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    k  = k+1;   
    resvec(k) = nr;
    % Compute normalized search direction. 
    if k==1
        p  = r;
        np = norm(p);
        p  = p/np; 
        q  = A*p;
        i  = 1; 
    else
        mk        = min(k-1, m);
        vim       = i+1:mk;
        ind       = [vim 1:i];
        P(:, i)   = p; 
        Q(:, i)   = q; 
        Qind = Q(:,ind);
        L(i, ind) = p'*Qind;
        v         = A*r;
        eP        = P(:,ind);  % Column sorting of P.
        lam       = tril(L(ind,ind)) \ (eP'*v);
        p         = r - eP*lam; 
        np        = norm(p);    
        p         = p/np; 
        q         = v - Qind*lam;  
        q         = q/np;           %q = A*p
        i         = i+1;
        if i == m+1
            i = 1;
        end
    end    
    %  Main step and residual update.
    pq = p'*q;
    a  = (nr^2)/(pq*np);         % step size
    x  = x + a*p;
    r  = r - a*q;
    nr = norm(r);
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end