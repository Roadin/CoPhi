%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%
% This code is a demo of the Maxide[1] method on synthetic data when the size
% of the matrix is 1000x1000, rank r=10, and r_A=r_B=20.
%
% [1] M. Xu, R. Jin, Z-H. Zhou. Speedup matrix completion with side
% information: application to multi-label data. In: NIPS, 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: Set Default Value and Randomization Seed
addpath('Maxide');
n=1000;
m=n;
r=10;
r_A=2*r;
r_B=2*r;
observed_rate=r*(2*n-r)/n/m;
lambda=1e-10;
max_iter=5000;

trial=1;
s = RandStream.create('mt19937ar','seed',trial);
%For MATLAB version 2013a or later, using
RandStream.setGlobalStream(s);
%For MATLAB version lower than 2013a, using
%RandStream.setDefaultStream(s);

disp('---------Matrix Information: ---------');
disp(['Matrix Size: ' num2str(n) ' x ' num2str(m) ])
disp(['Matrix Rank: ' num2str(r)]);
disp(['r_A = ' num2str(r_A)]);
disp(['r_B = ' num2str(r_B)]);
disp(['Observed Rate: ' num2str(observed_rate)]);
disp(['The regularization parameter is ' num2str(lambda)]);

%% STEP 1: Generate ground truth matrix M,A,B and Omega
disp('---------Matrix Generation: ---------');
L = randn(n,max(r_A,r_B));
R = randn(m,max(r_A,r_B));  

[A,~,B]=svds(L*R',max(r_A,r_B));
if n>5000 || r>10
    warning('For large matrix, the PROPACK (http://soi.stanford.edu/~rmunk/PROPACK/) is recommended for doing SVD:[A,~,B]=lansvd(L*R'',max(r_A,r_B),''L'')');
end
%addpath('PROPACK');
%[A,~,B]=lansvd(L*R',max(r_A,r_B),'L');


A=A(:,1:r_A);
B=B(:,1:r_B);
 
rand_U=randn(r_A,r);
U=A*rand_U;
[~,U_S,U_V]=svd(U,'econ');
multi_U=rand_U/(U_V')/U_S;
 
rand_V=randn(r_B,r);
V=B*rand_V;
[~,V_S,V_V]=svd(V,'econ');
multi_V=rand_V/(V_V')/V_S;
 
S=diag(randn(1,r));
 
trueZ=multi_U*S*multi_V';
trueZ=trueZ*(10^4);
M=A*trueZ*B';

Omega_linear = randsample(m*n, floor(m*n*observed_rate));
M_Omega=zeros(n,m);
M_Omega(Omega_linear)=M(Omega_linear);
M_Omega=sparse(M_Omega);
%% STEP 2: Recover M from Omega_linea, A and B
disp('---------Matrix Recovery: ---------');
[M_recover_side,telapsed_side]=Maxide(M_Omega,Omega_linear,A,B,lambda,max_iter);
l2loss_side= norm(M_recover_side-M,'fro')/norm(M,'fro');

disp('---------Recovery Result: ---------');
disp(['The relative error is: ' num2str(l2loss_side) ', ']);
disp(['The training time is ' num2str(telapsed_side) 'seconds.'])

