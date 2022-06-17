%% 2021-06-01
%% Na Huang

%  SWI   Sliding window implementation with pre-allocated memory
% 
%   [Z] = SWI(A, B) attempts to find a solution Z to the system of linear equations
%   AZ=B.  The N-by-N coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector B must have length N.
%
%   [Z] = SWI(A, B, MK) specifies the number of the sliding window. If MK is []
%   then SWI2 uses the default, N.
%
%   [Z] = SWI(A, B, MK, TOL) specifies the tolerance of the method. If TOL is []
%   then SWI2 uses the default, 1e-6.
%
%   [Z] = SWI(A, B, MK, TOL, Z)  specifies the initial guess.  If Z is [] then 
%   SWI2 uses the default, an all zero vector.
%
%   [Z] = SWI(A, B, MK, TOL, Z, MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then SWI2 uses the default, 10000.
%
%   [Z, K]=SWI(A, B, ...) returns the iteration number at which Z
%   was computed: 1 <= ITER <= MAXIT.
%
%   [Z, K, RES]=SWI(A, B, ...) also returns the last relative
%   residual norm NORM(B-AZ)/NORM(B).
%
%   [Z, K, Res, RESVEC]=SWI(A, B, ...) also returns a vector of estimates of the 
%   residual norms at each iteration, including NORM(B-AZ).
%------------------------------------------------------------------

function [x, k, res, resvec] = swi(A, b, mk, tol, x0, maxit)

n  = length(b);    
nb = norm(b);


%------------------------------------------------------------------
%     Retrieve input arguments.
%-----

if (nargin < 2)
   error('MATLAB:NotEnoughInputs', 'Not enough input arguments.');
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


L      = zeros(mk,mk);
P      = zeros(n,mk);
Q      = zeros(n,mk); 
resvec = zeros(maxit,1);



k   = 0;
r   = b - A*x;  
nr  = norm(r);
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
        mkk       = min(k-1, mk);
        vim       = i+1:mkk;
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
        if i == mk+1
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
