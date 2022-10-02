function [x, k, res, resvec] = diom(A, b, m, tol, x0, maxit)

% DIOM Direct Incomplete Orthogonalization Method
% Correspond to Algorithms 6.6 and 6.8 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"

%   [x] = diom(A, b) attempts to find a solution x to the system of linear equations
%   Ax=b. The n-by-n coefficient matrix A must be positive definite but need
%   not be symmetric. The right hand side column vector b must have length n.
%
%   [x] = diom(A, b, m) specifies the number of the sliding window. If m is [] 
%   then diom uses the default, n.
%
%   [x] = diom(A, b, m, tol) specifies the tolerance of the method. If tol is []
%   then diom uses the default, 1e-6.
%
%   [x] = diom(A, b, m, tol, x0)  specifies the initial guess. If x0 is [] 
%   then diom uses the default, an all zero vector.
%
%   [x] = diom(A, b, m, tol, x0, maxit) specifies the maximum number of iterations.
%   If maxit is [] then diom uses the default, 10000.
%
%   [x, k] = diom(A, b, ...) returns the iteration number at which x
%   was computed: 1 <= k <= maxit.
%
%   [x, k, res] = diom(A, b, ...) also returns the last relative
%   residual norm norm(b-Ax)/norm(b).
%
%   [x, k, res, resvec] = diom(A, b, ...) also returns a vector of estimates of the 
%   residual norms at each iteration, including norm(b-Ax).
%------------------------------------------------------------------
% 5 Jun 2021
% Na Huang

n  = length(b);


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
nb = norm(b);
x  = zeros(n,1);

mm1  = m-1;

resvec = zeros(maxit, 1);
%V      = zeros(n, m);
%P      = zeros(n, mm1);
U      = zeros(m,1);
L      = zeros(mm1,1);
h      = zeros(m+2,1); %% store the last column of Hk, h_{k+1,k}, and h_{k,k-1}


r      = b;
nr     = nb;
beta   = nr;
zeta   = beta;
V(:,1) = r/beta;




tol   = tol*nb;

k     = 0;

vm   = 2:m;
vmm1 = 2:mm1;
vmm2 = vmm1-1;
vmm3 = 1:mm1;

%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    k = k+1;
    resvec(k) = nr;
    if k==1
        V1 = V(:,1);
        w = A*V1;
        h(1) = w'*V1;
        w = w - h(1)*V1;
        h(2) = norm(w);
        V(:,2) = w/h(2);
        U(1) = h(1);
        p = V1/U(1);
        P(:,1) = p;
        x = x + zeta*p;
    else if k < m
            %% Incomplete Orthogonalization Process
            w = A*V(:,k);
            h(k+2) = h(k);
            for j = 1:k
                Vj   = V(:,j);
                h(j) = w'*Vj;
                w    = w - h(j)*Vj;
            end
            h(k+1) = norm(w);
            V(:,k+1) = w/h(k+1);
            %% Update the LU factorization of Hm
            L(k-1) = h(k+2)/U(k-1);
            U(1) = h(1);
            for i = 2:k-1
                U(i) = h(i) - L(i-1)*U(i-1);     % Compute U(:,k) = L\H(:,k)
            end
            U(k) = h(k) - L(k-1)*U(k-1);
            zeta = -L(k-1)*zeta;
            p = (V(:,k)-P(:,1:k-1)*U(1:k-1))/U(k);
            P(:,k) = p;
            x = x + zeta*p;
        else
            if k == m
                h(m+1) = h(m);
                U(m) = U(mm1);
            end
            %% Incomplete Orthogonalization Process
            w = A*V(:,m);
            h(m+2) = h(m+1);
            for j = 1:m
                Vj   = V(:,j);
                h(j) = w'*Vj;
                w    = w - h(j)*Vj;
            end
            h(j+1) = norm(w);
            V = [V(:,vm), w/h(j+1)];    % store V=[v_{k-m+2},...,v_{k},v_{k+1}]
            %% Update the LU factorization of Hm
            L(mm1) = h(m+2)/U(m);
            U(1) = h(1);
            for i = vmm1
                U(i) = h(i) - L(i-1)*U(i-1);     % Compute U(:,m) = L\H(:,m)
            end
            U(m) = h(m) - L(mm1)*U(mm1);
            zeta = -L(mm1)*zeta;
            p = (V(:,mm1)-P*U(vmm3))/U(m);
            x = x + zeta*p;
            P = [P(:,vmm1), p];
            L(vmm2) = L(vmm1);
        end
    end
    r  = b - A*x;
    nr = norm(r);
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end
