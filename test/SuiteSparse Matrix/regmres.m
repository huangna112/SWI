function [x, k, res, resvec] = regmres(A, b, restart, tol, x0, maxit)

% resgmres Restarted GMRES
% Correspond to Algorithms 6.6 and 6.11 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"

%   [x] = resgmres(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b. The n-by-n coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector b must have length n.
%
%   [x] = resgmres(A, b, restart) specifies the restarted number. If restart is []
%   then rescg uses the default, n.
%
%   [x] = resgmres(A, b, restart, tol) specifies the tolerance of the method. If tol is []
%   then resgmres uses the default, 1e-6.
%
%   [x] = resgmres(A, b, restart, tol, x0)  specifies the initial guess. If x0 is []
%   then resgmres uses the default, an all zero vector.
%
%   [x] = resgmres(A, b, restart, tol, x0, maxit) specifies the maximum number of iterations.
%   If maxit is [] then resgmres uses the default, 10000.
%
%   [x, k] = resgmres(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = resgmres(A, b, ...) also returns the last relative
%   residual norm norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = resgmres(A, b, ...) also returns a vector of estimates of the
%   residual norms at each iteration, including norm(b-Ax).
%------------------------------------------------------------------
% 2022-09-16
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

%% sovle Ad = A(x-x0) = b-Ax0 = r0 first
if gotx0 ~= 0
    b = b - A*x0;
end
nb = norm(b);
x  = zeros(n,1);


resvec = zeros(maxit, 1);
s      = zeros(restart,1);
c      = zeros(restart,1);
h      = zeros(restart+1,1); %% store the last column of Hk and h_{k+1,k}.




r      = b;
nr     = nb;
gamma0 = nr;
V(:,1) = r/gamma0;
tol    = tol*nb;
k      = 0;



%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    for inner = 1:restart
        k = k+1;
        resvec(k) = nr;
        %% Orthogonalization Process
        w = A*V(:,inner);
        for j = 1:inner
            Vj = V(:,j);
            h(j) = w'*Vj;
            w = w - h(j)*Vj;
        end
        h(inner+1) = norm(w);
        V(:,inner+1) = w/h(inner+1);
        %% Update the QR factorization of barHm
        if inner > 1
            for i = 1:inner-1
                temp = c(i)*h(i)+s(i)*h(i+1);
                h(i+1) = -s(i)*h(i)+c(i)*h(i+1);
                h(i) = temp;
            end
        end
        sqh = sqrt(h(inner)^2+h(inner+1)^2);
        s(inner) = h(inner+1)/sqh;
        c(inner) = h(inner)/sqh;
        gamma = -s(inner)*gamma0;
        gamma0 = c(inner)*gamma0;
        if inner==1
            p = V(:,inner)/sqh;
        else
            p = (V(:,inner)-P(:,1:inner-1)*h(1:inner-1))/sqh;
        end
        P(:,inner) = p;
        x = x + gamma0*p;
        nr  = abs(gamma);
        if (nr < tol)
            break;
        end
        gamma0 = gamma;
    end
    r      = b - A*x;
    gamma0 = nr;
    V(:,1) = r/gamma0;
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end