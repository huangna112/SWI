function [x, k, res, resvec] = diom(A, b, mk, tol, x0, maxit)

% Direct Incomplete Orthogonalization Method (DIOM)
% Correspond to Algorithms 6.6 and 6.8 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"
% 5 Jun 2021
% Na Huang

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

mkm1  = mk-1;

resvec = zeros(maxit, 1);
%V      = zeros(n, mk);
%P      = zeros(n, mkm1);
U      = zeros(mk,1);
L      = zeros(mkm1,1);
h      = zeros(mk+2,1); %% store the last column of Hk, h_{k+1,k}, and h_{k,k-1}


r      = b - A*x;
nr     = norm(r);
beta   = nr;
zeta   = beta;
V(:,1) = r/beta;




tol   = tol*nb;

k     = 0;

vmk   = 2:mk;
vmkm1 = 2:mkm1;
vmkm2 = vmkm1-1;
vmkm3 = 1:mkm1;

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
    else if k < mk
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
            if k == mk
                h(mk+1) = h(mk);
                U(mk) = U(mkm1);
            end
            %% Incomplete Orthogonalization Process
            w = A*V(:,mk);
            h(mk+2) = h(mk+1);
            for j = 1:mk
                Vj   = V(:,j);
                h(j) = w'*Vj;
                w    = w - h(j)*Vj;
            end
            h(j+1) = norm(w);
            V = [V(:,vmk), w/h(j+1)];    % store V=[v_{k-mk+2},...,v_{k},v_{k+1}]
            %% Update the LU factorization of Hm
            L(mkm1) = h(mk+2)/U(mk);
            U(1) = h(1);
            for i = vmkm1
                U(i) = h(i) - L(i-1)*U(i-1);     % Compute U(:,mk) = L\H(:,mk)
            end
            U(mk) = h(mk) - L(mkm1)*U(mkm1);
            zeta = -L(mkm1)*zeta;
            p = (V(:,mkm1)-P*U(vmkm3))/U(mk);
            x = x + zeta*p;
            P = [P(:,vmkm1), p];
            L(vmkm2) = L(vmkm1);
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