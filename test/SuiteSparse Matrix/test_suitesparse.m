fid =fopen('testproblem.txt','r');  % read file that contains test problems
sfil='TestResults';

leg={'BICGSTAB','GMRES','DQGMRES(5)','DQGMRES(10)','DQGMRES(100)','FOM','DIOM(2)','DIOM(5)','DIOM(10)','DIOM(100)','SCG','SWI(2)','SWI(5)','SWI(10)','SWI(100)'...
...'RGMRES', 'RFOM', 'RSCG', 'DSWI'};

% Set input parameters:
tol = 1.e-6;

% Initialize storage for output information

TestResult = cell(2,7);    %% name; size;
%% title
TestResult{1, 1} = 'name';
TestResult{1, 2} = 'size';
i = 4;
TestResult{1, i} = leg(1);
i = i+3;
TestResult{1, i} = leg(2);
i = i+3;
TestResult{1, i} = leg(3);
i = i+3;
TestResult{1, i} = leg(4);
i = i+3;
TestResult{1, i} = leg(5);
i = i+3;
TestResult{1, i} = leg(6);
i = i+3;
TestResult{1, i} = leg(7);
i = i+3;
TestResult{1, i} = leg(8);
i = i+3;
TestResult{1, i} = leg(9);
i = i+3;
TestResult{1, i} = leg(10);
i = i+3;
TestResult{1, i} = leg(11);
i = i+3;
TestResult{1, i} = leg(12);
i = i+3;
TestResult{1, i} = leg(13);
i = i+3;
TestResult{1, i} = leg(14);
i = i+3;
TestResult{1, i} = leg(15);
i = i+3;
TestResult{1, i} = leg(16);
i = i+3;
TestResult{1, i} = leg(17);
i = i+3;
TestResult{1, i} = leg(18);
i = i+3;
TestResult{1, i} = leg(19);

maxit = 10000;

p=2;
while ~feof(fid)
    name = fgetl(fid);
    load(name)
    %% Output the name
    TestResult{p, 1} = name;
    %% produce the test linear system
    A = Problem.A;
    n = size(A,1);
    b = A*ones(n,1);
    %% Output the size
    TestResult{p, 2} = n;
    
    
    
    %% BICGSTAB
    i = 3;
    tic
    [~, ~, TestResult{p, i+2}, TestResult{p, i}]=bicgstab(A, b,tol,maxit);
    TestResult{p, i+1} = toc;
    
    
    %% GMRES
    i = i+3;
    tic
    [~,TestResult{p, i},TestResult{p, i+2}] = regmres(A, b, [], tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = dqgmres(A, b, 5, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = dqgmres(A, b, 10, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = dqgmres(A, b, 100, tol);
    TestResult{p, i+1} = toc;
    
    
    %% FOM
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = refom(A, b, [], tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = diom(A, b, 2, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = diom(A, b, 5, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = diom(A, b, 10, tol);
    TestResult{p, i+1} = toc;
    
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}] = diom(A, b, 100, tol);
    TestResult{p, i+1} = toc;
    
    %% SCG
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=rescg(A, b, [], tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=swi(A, b, 2, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=swi(A, b, 5, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=swi(A, b, 10, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=swi(A, b, 100, tol);
    TestResult{p, i+1} = toc;
    
    %% Restart
    restart = 20;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=regmres(A, b, restart, tol);
    TestResult{p, i+1} = toc;
    
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=refom(A, b, restart, tol);
    TestResult{p, i+1} = toc;
        
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=rescg(A, b, restart, tol);
    TestResult{p, i+1} = toc;
    
    %% DSWI
    i = i+3;
    tic
    [~, TestResult{p, i}, TestResult{p, i+2}]=dynswiwp(A, b, restart, tol);
    TestResult{p, i+1} = toc;
    
    p = p+1
    save(sfil);
end

