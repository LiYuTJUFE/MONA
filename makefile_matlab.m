clear
cpath = pwd;
ipath = ['-I' cpath];
% mex -setup C
% int in C should be replaced by long long int, which is 8 bytes.
mex('-v','-R2017b','-Dint="long long int"',ipath,'gcge_matlab.c',...
    'ops.c','ops_multi_vec.c',...
    'ops_orth.c','ops_multi_grid.c',...
    'ops_lin_sol.c','ops_eig_sol_gcg.c','ops_eig_sol_pas.c',...    
    'app_lapack.c','app_ccs.c','app_slepc.c','app_pas.c')
%%
clear
format long
n=200; h = 1/(n+1);
A=sparse([1:(n-1),1:n,2:n],[2:n,1:n,1:(n-1)],[-1*ones(1,n-1),2*ones(1,n),-1*ones(1,n-1)],n,n);
A=A/h;%full(A)
B=sparse([1:n],[1:n],[1*ones(1,n)],n,n);
B=B*h;%full(B)
nev=10; multiMax=2; gapMin=0.01; block_size=3; tol=[1e-2;1e-8]; numIterMax=100;
[eval,evec,nevConv] = gcge_matlab(A,B,nev,multiMax,gapMin,block_size,tol,numIterMax);
meval=eigs(A,B,nev,'smallestabs');
eval(1:nev)-meval(1:nev)
