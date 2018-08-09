function [G_new, E_new] = TCSIM(S_1, S_2, S_3, T, lambda_T, lambda_C, lambda_N, maxIter, eps, zeta)
%Input: X -- User matrix
%       Y -- Movie matrix
%       Z -- Time matrix
%     
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
 nonzeros = F>0;
 Omega = double(nonzeros); 
 beta_max = 100;
 rho = 1.01;
 beta_0 = 4*1e-5;

 dump_percent = 0.1;
dim = size( F ); % dimension of the matrix we want to complete

%% p is features of matrix X; q is features of matrix Y
d_1 = size(X, 1);
d_2 = size(Y, 1);
d_3 = size(Z, 1);
G = zeros(d_1, d_2, d_3);

%% First, we use SVT  to find initial matrix
temp = reshape( F, m * n *p, 1);
[Tau, col] = find( ~(temp == 0));
data = temp( Tau );
tau = 5*sqrt(m * n * p); 
delta = 1;
maxit = floor( 500/dump_percent);
tol = 1e-1;
dim = size( F );
F_mode_1 = Unfold(F, dim, 1);
[U, S, V, numiter] = SVT(size(F_mode_1), Tau, data, tau, delta, maxit, tol);
E =  U*S*V';
E(F~=0) = F(F~=0);
E = Fold(E, dim, 1);


multiplier_1 = zeros(size(E));
multiplier_2 = zeros(size(E));

%% ADM Iteration 
obj = 0;
A = kron(Y', X'); %kronecker tensor product
A = kron(Z', A);
beta = beta_0;
G = Kron_G(A, Unfold(E, dim, 1), lambda_G, dim(1), dim(2));
G = Fold(G, [d_1, d_2, d_3], 1);

for iter = 1 : maxIter
    %update C
    C_new = closed_form_C(beta, E, G, X, Y, Z, multiplier_2);
    % update zeta
    % update G
    G_new = soft_one_norm_G(A, E, G, C_new, lambda_G, multiplier_2, zeta, beta);
    G_new = Fold(G_new, [d_1, d_2, d_3], 1);
    %  update E
    E_new = singular_value_th_E(F, E, G_new, C_new, Omega, X, Y, Z, lambda_E, multiplier_1, multiplier_2, zeta, beta);
    % update M_1, M_2
    [multiplier_1_new, multiplier_2_new] = multiplier_update(multiplier_1, multiplier_2, E_new, F, G_new, X, Y, C_new, Omega, beta);
    % update beta
    beta_set = [beta_max, rho*beta];
    beta_new = min(beta_set);
    %objective function value
    nulear_norm_E = 0;
    for dim_i = 1 : 3
        e = svd(Unfold(E_new, dim, 1));
        nulear_norm_E = sum(e) + nulear_norm_E;
    end
        
    obj_new = lambda_G * norm(G_new(:), 1) + lambda_E * (nulear_norm_E) + 0.5 * norm(Unfold(C_new, dim, 1), 'fro')^2;
    
    fprintf('TCSIM: iterations = %d   objective funciton = %f\n', iter, obj_new);
    
    err_obj_per = abs( obj_new - obj)/abs(obj);
    
    if err_obj_per < eps
        break;
    end
    
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


end