format short e
tol=1.e-6;
load A.mat
load b.mat
A = Asupg;
b = fsupg;



fid = 1;

maxit = 10000;

%% BICGSTAB
tic
[z_big, flag_big, res_big, k_big, resvec_big]=bicgstab(A, b,tol,maxit);
t_big = toc;
fprintf(fid, 'BICGSTAB,             %4i    %4.4f    %11.4e  \n',k_big*2,t_big,res_big); fprintf('\n');

%% GMRES
tic
[z_gmres,k_gmres,res_gmres, resvec_gmres] = dqgmres(A, b, [], tol);
t_gmres = toc;
fprintf(fid, 'GMRES,                %4i     %4.4f     %14.4e \n', k_gmres, t_gmres, res_gmres); fprintf('\n');

% tic
% [z_dqgmres, k_dqgmres, res_dqgmres, resvec_dqgmres] = dqgmres(A, b, 100, tol);
% t_dqgmres = toc;
% fprintf(fid, 'DQGMRES(100),          %4i     %4.4f     %14.4e \n', k_dqgmres, t_dqgmres, res_dqgmres); fprintf('\n');


%% FOM
tic
[z_fom, k_fom, res_fom, resvec_fom] = diom(A, b, [], tol);
t_fom = toc;
fprintf(fid, 'FOM,                  %4i    %4.4f    %11.4e  \n',k_fom,t_fom,res_fom); fprintf('\n');

tic
[z_diom1, k_diom1, res_diom1, resvec_diom1] = diom(A, b, 2, tol);
t_diom1 = toc;
fprintf(fid, 'DIOM(2),             %4i    %4.4f    %11.4e  \n',k_diom1,t_diom1,res_diom1); fprintf('\n');

tic
[z_diom2, k_diom2, res_diom2, resvec_diom2] = diom(A, b, 5, tol);
t_diom2 = toc;
fprintf(fid, 'DIOM(5),             %4i    %4.4f    %11.4e  \n',k_diom2,t_diom2,res_diom2); fprintf('\n');

tic
[z_diom3, k_diom3, res_diom3, resvec_diom3] = diom(A, b, 10, tol);
t_diom3 = toc;
fprintf(fid, 'DIOM(10),             %4i    %4.4f    %11.4e  \n',k_diom3,t_diom3,res_diom3); fprintf('\n');


tic
[z_scg, k_scg, res_scg, resvec_scg]=swiwp(A, b, [], tol);
t_scg = toc;
fprintf(fid, 'SCG,                  %4i    %4.4f    %11.4e  \n',k_scg,t_scg,res_scg); fprintf('\n');

tic
[z_swi1, k_swi1, res_swi1, resvec_swi1]=swi(A, b, 2, tol);
t_swi1 = toc;
fprintf(fid, 'SWI(2),              %4i    %4.4f    %11.4e  \n',k_swi1,t_swi1,res_swi1); fprintf('\n');

tic
[z_swi2, k_swi2, res_swi2, resvec_swi2]=swi(A, b, 5, tol);
t_swi2 = toc;
fprintf(fid, 'SWI(5),              %4i    %4.4f    %11.4e  \n',k_swi2,t_swi2,res_swi2); fprintf('\n');

tic
[z_swi3, k_swi3, res_swi3, resvec_swi3]=swi(A, b, 10, tol);
t_swi3 = toc;
fprintf(fid, 'SWI(10),              %4i    %4.4f    %11.4e  \n',k_swi3,t_swi3,res_swi3); fprintf('\n');


