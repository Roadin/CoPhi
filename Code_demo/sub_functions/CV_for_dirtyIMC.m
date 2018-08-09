function [training_error, Elapsed_time] = CV_for_dirtyIMC(dirtyIMC_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix)
    [var_set, var_size ] = generate_variable_permutation(dirtyIMC_params.lambda_set, dirtyIMC_params.lambda_1_set, dirtyIMC_params.iteration, [-1], [-1]);
    relative_error_matrix_for_each_method = ones(var_size,1);

    for i=1:var_size
        lambda = var_set(i,1);
        lambda_1 = var_set(i,2);
        maxiter = var_set(i,3);
        [UU, SS, VV, U, S, V] = dirtyIMC(training_Data_inside, X', Y', lambda, lambda_1, maxiter,0,1e-6);
        M = UU*SS*VV'; 
        N = U*S*V'; 
        E_new = X'*M*Y + N;
        [relative_error_matrix_for_each_method(i,1),~] = calculate_error(E_new, training_Data, dump_position);
    end
    best_relative_error = min(relative_error_matrix_for_each_method);
    validating_var_set = var_set(find(relative_error_matrix_for_each_method == best_relative_error(1)),:);

    k = validating_var_set(1, 1);
    lambda = validating_var_set(1, 2);
    maxiter = validating_var_set(1, 3);
%     W0 = randn(d1, k);
%     H0 = randn(d2, k);
    tic;
    [UU, SS, VV, U, S, V] = dirtyIMC(training_Data, X', Y', lambda, lambda_1, maxiter,0,1e-6);
    M = UU*SS*VV'; 
    N = U*S*V'; 
    E_new = X'*M*Y + N;
    Elapsed_time = toc;

    [training_error, ~] = calculate_error(E_new, testing_Data, dump_matrix);    
end