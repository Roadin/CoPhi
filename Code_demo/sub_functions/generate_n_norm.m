function Y = generate_n_norm(E, dim, n)
    Y = zeros(dim, dim);
    
    for i = 1:dim
        for j = 1:dim
            a = E(:,i);
            b = E(:,j);
            if n == -1
                Y(i,j) = a'*b/norm(a)/norm(b);
            else
                Y(i,j) = norm(a-b, n);
            end
        end
    end
    if i ~= -1
        Y = Y ./ max(Y(:));
        Y = 1 - Y;
    end
    clear i j a b current_norm
end