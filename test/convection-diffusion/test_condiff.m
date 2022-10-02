%% Test of the linear systems from the convection-diffusion equations

%------------------------------------------------------------------
%     Initialization 
%-----
format short e
tol   = 1.e-6;
maxit = 10000;
fid   = 1;
% input test data
load A.mat
load b.mat
A = Asupg;
b = fsupg;

%% BICGSTAB
tic
[z_big, flag_big, res_big, k_big, resvec_big] = bicgstab(A, b,tol,maxit);
t_big = toc;
fprintf(fid, 'BICGSTAB,         %4i    %4.4f     %11.2e  \n',k_big*2+1,t_big,res_big); fprintf('\n');

%% GMRES
tic
[z_gmres,k_gmres,res_gmres, resvec_gmres] = regmres(A, b, [], tol);
t_gmres = toc;
fprintf(fid, 'GMRES,            %4i    %4.4f    %11.2e  \n', k_gmres, t_gmres, res_gmres); fprintf('\n');

%% DQGMRES
tic
[z_dqgmres,k_dqgmres,res_dqgmres, resvec_dqgmres] = dqgmres(A, b, 36, tol);
t_dqgmres = toc;
fprintf(fid, 'DQGMRES(36),      %4i    %4.4f    %11.2e  \n', k_dqgmres, t_dqgmres, res_dqgmres); fprintf('\n');

%% FOM
tic
[z_fom, k_fom, res_fom, resvec_fom] = refom(A, b, [], tol);
t_fom = toc;
fprintf(fid, 'FOM,              %4i    %4.4f    %11.2e  \n',k_fom,t_fom,res_fom); fprintf('\n');


%% DIOM(mk)
tic
[z_diom1, k_diom1, res_diom1, resvec_diom1] = diom(A, b, 2, tol);
t_diom1 = toc;
fprintf(fid, 'DIOM(2),          %4i    %4.4f    %11.2e  \n',k_diom1,t_diom1,res_diom1); fprintf('\n');

tic
[z_diom2, k_diom2, res_diom2, resvec_diom2] = diom(A, b, 5, tol);
t_diom2 = toc;
fprintf(fid, 'DIOM(5),          %4i    %4.4f    %11.2e  \n',k_diom2,t_diom2,res_diom2); fprintf('\n');

tic
[z_diom3, k_diom3, res_diom3, resvec_diom3] = diom(A, b, 10, tol);
t_diom3 = toc;
fprintf(fid, 'DIOM(10),          %4i    %4.4f    %11.2e  \n',k_diom3,t_diom3,res_diom3); fprintf('\n');

%% SCG
tic
[z_scg, k_scg, res_scg, resvec_scg] = rescg(A, b, [], tol);
t_scg = toc;
fprintf(fid, 'SCG,               %4i    %4.4f    %11.2e  \n',k_scg,t_scg,res_scg); fprintf('\n');


%% SWI(mk)
tic
[z_swi1, k_swi1, res_swi1, resvec_swi1] = swi(A, b, 2, tol);
t_swi1 = toc;
fprintf(fid, 'SWI(2),            %4i    %4.4f    %11.2e  \n',k_swi1,t_swi1,res_swi1); fprintf('\n');

tic
[z_swi2, k_swi2, res_swi2, resvec_swi2] = swi(A, b, 5, tol);
t_swi2 = toc;
fprintf(fid, 'SWI(5),            %4i    %4.4f    %11.2e  \n',k_swi2,t_swi2,res_swi2); fprintf('\n');

tic
[z_swi3, k_swi3, res_swi3, resvec_swi3] = swi(A, b, 10, tol);
t_swi3 = toc;
fprintf(fid, 'SWI(10),           %4i    %4.4f    %11.2e  \n',k_swi3,t_swi3,res_swi3); fprintf('\n');


%% Restart

Restart = 20;

%% RGMRES
tic
[z_rgmres,k_rgmres,res_rgmres, resvec_rgmres] = regmres(A, b, Restart, tol);
t_rgmres = toc;
fprintf(fid, 'RGMRES,           %4i    %4.4f    %11.2e  \n', k_rgmres, t_rgmres, res_rgmres); fprintf('\n');


%% RFOM
tic
[z_rfom, k_rfom, res_rfom, resvec_rfom] = refom(A, b, Restart, tol);
t_rfom = toc;
fprintf(fid, 'RFOM,             %4i    %4.4f    %11.2e  \n',k_rfom,t_rfom,res_rfom); fprintf('\n');


%% RSCG
tic
[z_rscg, k_rscg, res_rscg, resvec_rscg] = rescg(A, b, Restart, tol);
t_rscg = toc;
fprintf(fid, 'RSCG,             %4i    %4.4f    %11.2e  \n',k_rscg,t_rscg,res_rscg); fprintf('\n');


%% DSWI
tic
[z_dyswi, k_dyswi, res_dyswi, resvec_dyswi] = dynswiwp(A, b, 2, tol);
t_dyswi = toc;
fprintf(fid, 'DSWI,             %4i    %4.4f    %11.2e  \n',k_dyswi,t_dyswi,res_dyswi); fprintf('\n');

%% draw res
k_big=1:(k_big*2+1);  k_gmres=1:k_gmres;  k_dqgmres=1:k_dqgmres;
k_fom=1:k_fom;  k_diom2=1:k_diom2; k_diom3=1:k_diom3; 
k_scg=1:k_scg;  k_swi1=1:k_swi1; k_swi2=1:k_swi2; k_swi3=1:k_swi3;
k_rgmres=1:k_rgmres; k_rfom=1:k_rfom; k_rscg=1:k_rscg; k_dyswi=1:k_dyswi;

figure(1)
semilogy(k_big,resvec_big,'r','Linewidth', 1.5); hold on
semilogy(k_gmres,resvec_gmres(k_gmres),'*-b','Linewidth', 1.0);
semilogy(k_dqgmres,resvec_dqgmres(k_dqgmres),'g','Linewidth', 1.5);
semilogy(k_fom,resvec_fom(k_fom),'m','Linewidth', 1.5);
semilogy(k_scg,resvec_scg(k_scg),'k','Linewidth', 1.5);
xlabel('Iter')
ylabel('Res')
legend('BICGSTAB','GMRES','DQGMRES(36)','FOM','SCG')
hold off


figure(2)
semilogy(k_diom2,resvec_diom2(k_diom2),'r','Linewidth', 1.5); hold on
semilogy(k_diom3,resvec_diom3(k_diom3),'b','Linewidth', 1.5);
semilogy(k_swi1,resvec_swi1(k_swi1),'g','Linewidth', 1.5);
semilogy(k_swi2,resvec_swi2(k_swi2),'m','Linewidth', 1.5);
semilogy(k_swi3,resvec_swi3(k_swi3),'k','Linewidth', 1.5);
xlabel('Iter')
ylabel('Res')
legend('DIOM(5)','DIOM(10)','SWI(2)','SWI(5)','SWI(10)')
hold off


figure(3)
semilogy(k_rgmres,resvec_rgmres(k_rgmres),'r','Linewidth', 1.5); hold on
semilogy(k_rfom,resvec_rfom(k_rfom),'b','Linewidth', 1.5);
semilogy(k_rscg,resvec_rscg(k_rscg),'g','Linewidth', 1.5);
semilogy(k_dyswi,resvec_dyswi(k_dyswi),'m','Linewidth', 1.5);
xlabel('Iter')
ylabel('Res')
legend('RGMRES','RFOM','RSCG','DSWI')
hold off
