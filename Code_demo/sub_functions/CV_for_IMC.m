function [training_error, Elapsed_time] = CV_for_IMC(IMC_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix)
    [var_set, var_size ] = generate_variable_permutation(IMC_params.k_set, IMC_params.lambda_set, IMC_params.iteration, [-1], [-1]);
    d1 = size( X,1);
    d2 = size( Y,1);
    relative_error_matrix_for_each_method = ones(var_size,1);

    for i=1:var_size
        k = var_set(i, 1);
        lambda = var_set(i, 2);
        maxiter = var_set(i, 3);
        W0 = randn(d1, k);
        H0 = randn(d2, k);

        [W, H] = IMC(training_Data_inside, X', Y', k, lambda, maxiter, W0, H0);
        E_new = X'*W'*H*Y;
        [relative_error_matrix_for_each_method(i,1),~] = calculate_error(E_new, training_Data, dump_position);    
    end
    best_relative_error = min(relative_error_matrix_for_each_method);
    validating_var_set = var_set(find(relative_error_matrix_for_each_method == best_relative_error),:);
    k = validating_var_set(1,1);
    lambda = validating_var_set(1, 2);
    maxiter = validating_var_set(1, 3);
    W0 = randn(d1, k);
    H0 = randn(d2, k);
    tic;
    [W, H] = IMC(training_Data, X', Y', k, lambda, maxiter, W0, H0);
    E_new = X'*W'*H*Y;
    Elapsed_time = toc;

    [training_error, ~] = calculate_error(E_new, testing_Data, dump_matrix);    
end