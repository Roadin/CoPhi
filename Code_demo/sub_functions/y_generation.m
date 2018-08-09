function [Y] = y_generation(Y_select, X_full_observe, F_full_observe)
    dim = size(F_full_observe,2);
    if strcmp(Y_select, 'Ya')    
    %% create Y_1
    Y = eye(dim,dim);
    
    elseif strcmp(Y_select, 'Yb')
        %% create diagonal matrix Y_2
        Y = eye(dim,dim);
        for i=1:11
            Y(i, i+11) = 1;
            Y(i+11, i) = 1;
        end
        clear i
    elseif strcmp(Y_select, 'Ycc')
        %% create similarity function Y_3c (cosine similarity)
        Y = generate_n_norm(F_full_observe, dim, -1);
    elseif strcmp(Y_select, 'Yci')
        %% create similarity function Y_3i (infinite norm)
        Y = generate_n_norm(F_full_observe, dim, Inf);
    elseif strcmp(Y_select, 'Y30n')
        %% create similarity function Y_30n (0 norm)
        % Y_30n = generate_n_norm(F_full_observe, 0);
    elseif strcmp(Y_select, 'Y31n')
        %% create similarity function Y_31n (1 norm)
        Y = generate_n_norm(F_full_observe, dim, 1);
    elseif strcmp(Y_select, 'Y32n')
        %% create similarity function Y_32n (2 norm)
        Y = generate_n_norm(F_full_observe, dim, 2);
    elseif strcmp(Y_select, 'Y33n')
        %% create similarity function Y_33n (2 norm)
        Y = generate_n_norm(F_full_observe, dim, 3);
    elseif strcmp(Y_select, 'Y4')
        %% create correlation function Y_4
        Y = corrcoef(F_full_observe);
        Y(isnan(Y)) = 0;
        Y(17,17) = 1;
    end
end