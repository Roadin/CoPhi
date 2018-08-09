function [G_new, E_new, Elapsed_time,iter] = AFA_ALDM_sto(X, Y, F, lambda_E, lambda_G, maxIter, eps,zeta, mult, beta_0)
%Input: X -- User matrix
%       Y -- Movie matrix
%       F -- Rating matrix ( m by n )
%       lambda_E, lambda_G -- the regularization parameter
%       maxIter -- maximum number of iterations
%       eps -- 
%       zeta -- control the step size: recommend 1e4 for small dataset,
%               1e10 to large dataset. When zeta is too small, the objective
%               function will not converge. When zeta is too large, the objective
%               function will stop before coverge. 

%Output: G_new -- the recovered G
%        E_new -- the recovered rating matrix E
%
%

%% parameters
nonzeros = F ~= 0;
Omega = double(nonzeros); 
beta_max = 100;
rho = 1.01;
 
dump_percent = 0.1;
[m, n] = size( F ); % dimension of the matrix we want to complete

%% p is features of matrix X; q is features of matrix Y
p = size(X, 1);
q = size(Y, 1);
G = zeros(p, q);

%% First, we use SVT  to find initial matrix
temp = reshape( F, m * n, 1 );
[Tau, col] = find( ~(temp == 0));
data = temp( Tau );
tau = 5*sqrt(m * n); 
delta = 1;
maxit = floor( 500/dump_percent);
tol = 1e-1;
[U, S, V, numiter] = SVT(size(F), Tau, data, tau, delta, maxit, tol);
E =  U*S*V';
E(F~=0) = F(F~=0);




multiplier_1 = zeros(size(E));
multiplier_2 = zeros(size(E));

%% ADM Iteration 
obj = 1e9;

[x_n1, x_n2] = size(X);

%% Improve kron by reducing size of X and F
[x_n1, x_n2] = size(X);
if x_n1 < x_n2
    A = kron(Y', X(1:x_n1, 1:x_n1)'); %kronecker tensor product
else
    A = kron(Y', X'); %kronecker tensor product    
end

beta = beta_0;
%     C_new = closed_form_C(beta, E, G, X, Y, multiplier_2);
%     G = soft_one_norm_G(A, E, G, C_new, lambda_G, multiplier_2, zeta, beta);
dim_G = size(G);

if x_n1 < x_n2
    G = Kron_G(A, E(1:x_n1, 1:n), lambda_G, p, q);
else
    G = Kron_G(A, E, lambda_G, p, q);
end
%     G = Initial_G_inv(A, E, dim_G);

% settings for parallel stochastic ladmm
max_worker_num = 8;
size_batch = ceil(size(A, 1)/mult);

tic;
for iter = 1 : maxIter
    % update C
    C_new = closed_form_C(beta, E, G, X, Y, multiplier_2);
    % update G
    G_new = soft_one_norm_stochastic_G(A, E, G, C_new, lambda_G, multiplier_2, zeta, beta, mult);
    %  update E
    E_new = singular_value_th_E(F, E, G_new, C_new, Omega, X, Y, lambda_E, multiplier_1, multiplier_2, zeta, beta);
    % update M_1, M_2
    [multiplier_1_new, multiplier_2_new] = multiplier_update(multiplier_1, multiplier_2, E_new, F, G_new, X, Y, C_new, Omega, beta);
    % update beta
    beta_set = [beta_max, rho*beta];
    beta_new = min(beta_set);
    % objective function value
    e = svd(E_new);
    nulear_norm_E = sum(e);
    
    
    obj_new = lambda_G * norm(G_new(:), 1) + lambda_E * nulear_norm_E + 0.5 * norm(C_new, 'fro')^2;
    
  
    
    err_obj_per = (obj-obj_new)/abs(obj_new);
    
    if (err_obj_per < eps || abs(obj-obj_new)<eps)
        break;
    end
     fprintf('ALDM: iterations = %d   objective funciton = %f\n', iter, obj_new); 
    if iter == maxIter
        fprintf('Not stop\n');
    end
    
    G = G_new;
    E = E_new;
    beta = beta_new;
    multiplier_1 = multiplier_1_new;
    multiplier_2 = multiplier_2_new;
    obj = obj_new;
end
Elapsed_time = toc;

end