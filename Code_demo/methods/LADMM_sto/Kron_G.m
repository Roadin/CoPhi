function new_G = Kron_G(A, E, lambda_G, p, q)
g = lasso(A, E(:), 'Lambda', 2 * lambda_G); %lasso fitted value
new_G = reshape(g, p, q);
end
