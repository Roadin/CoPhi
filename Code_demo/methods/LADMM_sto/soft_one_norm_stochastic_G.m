function G_new = soft_one_norm_stochastic_G(A, E, G, C_new, lambda_G, multiplier_2, tau, beta,mult)
%% Update G via linearization of the subproblem

%
% B = E + multiplier_2/beta;
% L = B - C_new;
 [G_row, G_col] = size(G);
% L = L(:);
% y = G(:);
% G = G(:);
% multiplier_2 = G(:);
% L_length = length(L);
% size_batch = 100;
% max_i = 100;
% threshold_G = lambda_G/beta;
% threshold_G = threshold_G *ones(length(G), 1);
% for i = 1 : max_i
%     picked_ind = randperm(L_length, size_batch);
%     y_new = 1/(beta+1/tau)*A(picked_ind, :)'*(L(picked_ind, :)-A(picked_ind, :)*y) + multiplier_2 + beta*G + y/tau;
%     G_new = soft_thresh(G-multiplier_2/beta, threshold_G);
%     multiplier_2 = multiplier_2 - beta*(y_new-G_new);
%     y = y_new;
%     G = G_new;
% end
% G_new = reshape(G, [G_row, G_col]);
% 
size_batch = ceil(size(A, 1)/mult);
picked_ind = randperm(size(A, 1), size_batch);
B = C_new(:) - E(:) - multiplier_2(:) / beta;
gradient_G = A(picked_ind, :)' * (A(picked_ind, :) * G(:) + B(picked_ind) );

%reshape matrix G to be vector and update
G_temp = G(:) - 1/(tau*beta) * gradient_G;
threshold_G =  lambda_G/(beta*tau) * ones(size(G(:)));
G_new = soft_thresh(G_temp, threshold_G);
G_new = reshape(G_new, [G_row, G_col]);
end