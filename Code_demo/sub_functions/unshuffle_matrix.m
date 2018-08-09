function [X_new,E_new] = unshuffled_matrix(X,E,shuffle_index)
    
    rows_E = size(E,1);
    X_new = X;
    E_new = E;
    
    X_new(:,shuffle_index)=X(:,[1:1:rows_E]);
    E_new(shuffle_index,:)=E([1:1:rows_E],:);
    clear E X shuffle_index
end

