function [error_rate, norm_value] = calculate_error(E_new, E_org, dump_matrix)
    elements = sum(sum(dump_matrix));
    dim = numel(E_new);
    error_matrix = abs(E_org-E_new).*dump_matrix;
    clear E_org E_filled dumped_matrix
    norm_value = norm(error_matrix, 'fro')/elements;
    error_matrix = double(error_matrix>=1);
    error_rate = sum(error_matrix(:))/elements;
end