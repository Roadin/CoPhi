function [X_new,E_new,shuffle_index] = shuffle_matrix(X,E)
    columns_X = size(X, 2);
    rows_E = size(E, 1);
    if columns_X ~= rows_E
        error('dimension not match')
    end
    clear columns_X
    
    shuffle_index = randperm(rows_E, rows_E);
    X_new = X;
    E_new = E;
    
    X_new(:,[1:1:rows_E])=X(:,shuffle_index);
    E_new([1:1:rows_E],:)=E(shuffle_index,:);
    clear E X i
end


