function [x, k, res, resvec] = dqgmres( A, b, mk, tol, x0, maxit)
% Direct Quasi-GMRES (DQGMRES)
% Correspond to Algorithms 6.6 and 6.13 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"
% 8 July 2021
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
s      = zeros(mk,1);
c      = zeros(mk,1);
h      = zeros(mk+1,1); %% store the last column of Hk and h_{k+1,k}.



r      = b - A*x;
nr     = norm(r);
gamma0 = nr;
V(:,1) = r/gamma0;
resvec(1) = nr;



tol   = tol*nb;

k     = 0;


vmk   = 2:mk;     %% [2 3 ... mk]
vmkm1 = 2:mkm1;   %% [2 3 ... mk-1]
vmkm3 = 1:mkm1;   %% [1 2 ... mk-1]

%------------------------------------------------------------------
%     Main iteration loop.
%-----
while (k < maxit) && (nr >= tol)
    k = k+1;
    if k==1
        V1 = V(:,1);
        w  = A*V1;
        h(1) = w'*V1;
        w = w - h(1)*V1;
        h(2) = norm(w);
        V(:,2) = w/h(2);
        sqh = sqrt(h(1)^2+h(2)^2);
        s(1) = h(2)/sqh;
        c(1) = h(1)/sqh;
        gamma = -s(1)*gamma0;
        gamma0 = c(1)*gamma0;
        p = V1/sqh;
        P(:,1) = p;
        x = x + gamma0*p;
        nr  = abs(gamma);
        resvec(k+1) = nr;
    else if k < mk
           %% Incomplete Orthogonalization Process
            w = A*V(:,k);
            for j = 1:k
                Vj = V(:,j);
                h(j) = w'*Vj;
                w = w - h(j)*Vj;
            end
            h(k+1) = norm(w);
            V(:,k+1) = w/h(k+1);
            %% Update the QR factorization of barHm
            for i = 1:k-1
                temp = c(i)*h(i)+s(i)*h(i+1);
                h(i+1) = -s(i)*h(i)+c(i)*h(i+1);
                h(i) = temp;
            end
            sqh = sqrt(h(k)^2+h(k+1)^2);
            s(k) = h(k+1)/sqh;
            c(k) = h(k)/sqh;
            gamma = -s(k)*gamma0;
            gamma0 = c(k)*gamma0;
            p = (V(:,k)-P(:,1:k-1)*h(1:k-1))/sqh;
            P(:,k) = p;
            x = x + gamma0*p;
            nr  = abs(gamma);
            resvec(k+1) = nr;
        else
            %% Incomplete Orthogonalization Process
            w = A*V(:,mk);
            for j = 1:mk
                Vj   = V(:,j);
                h(j) = w'*Vj;
                w    = w - h(j)*Vj;
            end
            h(mk+1) = norm(w);
            V = [V(:,vmk), w/h(mk+1)];    % store V=[v_{k-mk+2},...,v_{k},v_{k+1}]
            %% Update the QR factorization of Hm
            for i = 1:mkm1
                temp = c(i)*h(i)+s(i)*h(i+1);
                h(i+1) = -s(i)*h(i)+c(i)*h(i+1);
                h(i) = temp;
            end
            sqh = sqrt(h(mk)^2+h(mk+1)^2);
            s(mk) = h(mk+1)/sqh;
            c(mk) = h(mk)/sqh;
            gamma = -s(mk)*gamma0;
            gamma0 = c(mk)*gamma0;
            p = (V(:,mkm1)-P*h(vmkm3))/sqh;
            x = x + gamma0*p;
            P = [P(:,vmkm1), p];
            s(vmkm3) = s(vmk);
            c(vmkm3) = c(vmk);
            r  = b - A*x;
            nr = norm(r);
            resvec(k+1) = nr;
        end
    end
    gamma0 = gamma;
end
k = k+1;
resvec(k) = nr;
res = nr/nb;
if gotx0 ~= 0
    x   = x + x0;
end