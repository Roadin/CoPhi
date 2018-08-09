function G_new = soft_one_norm_G(A, E, G, C_new, lambda_G, multiplier_2, tau, beta)
%% Update G via linearization of the subproblem

%
gradient_G = A' * (A * G(:) + C_new(:) - E(:) - multiplier_2(:) / beta );

%reshape matrix G to be vector and update
G_temp = G(:) - 1/(tau*beta) * gradient_G;
threshold_G =  lambda_G/(beta*tau) * ones(size(G(:)));
G_new = soft_thresh(G_temp, threshold_G);

[G_row, G_col] = size(G);
G_new = reshape(G_new, [G_row, G_col]);
end