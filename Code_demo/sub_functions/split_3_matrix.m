function [train_matrix, valid_matrix, test_matrix] = split_3_matrix(E)
    rows = size(E,1);
    train_matrix = E(1:floor(rows/3), :);
    valid_matrix = E(floor(rows/3)+1:floor(2*rows/3), :);
    test_matrix = E(floor(2*rows/3)+1:rows, :);
end
