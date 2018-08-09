function new_G = Initial_g_inv(A, T, dim_G)
    G = inv(A'*A + 0.0001 .* eye(size(A'*A)))*A'*T(:);
    new_G = reshape(G, dim_G);
end
