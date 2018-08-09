function [training_error, Elapsed_time,G_new] = CV_for_LADMM(LADMM, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix)
    [var_set, var_size ] = generate_variable_permutation(LADMM.lambda_E_set, LADMM.lambda_G_set,LADMM.zeta_set,LADMM.mult_set,LADMM.beta_0_set);
    maxIter = 3000;
    eps = 10e-6;
    relative_error_matrix_for_each_method = ones(var_size,1);

    for i=1:var_size
        try
            [G_new, E_new,Elapsed_time] = AFA_ALDM(X, Y, training_Data_inside, var_set(i,1), var_set(i,2), maxIter, eps, var_set(i,3), var_set(i,5));
            E_XGY = X' * G_new * Y;
            
        [relative_error1, norm_value1] = calculate_error(E_new, training_Data, dump_position);
        [relative_error2, norm_value2] = calculate_error(E_XGY, training_Data, dump_position);
        relative_error_matrix_for_each_method(i,1) = min(relative_error1, relative_error2);
        catch
            warning('Error in LADMM');
            relative_error_matrix_for_each_method(i,1) = 1024;            
        end
        csvwrite( strcat(pwd, '//../data_journal/current_error_matrix.csv'), relative_error_matrix_for_each_method);        
    end
    best_relative_error = min(relative_error_matrix_for_each_method);
    validating_var_set = var_set(find(relative_error_matrix_for_each_method == best_relative_error),:);

    [G_new, E_new,Elapsed_time] = AFA_ALDM(X, Y, training_Data, validating_var_set(1,1), validating_var_set(1,2), maxIter, eps, validating_var_set(1,3), validating_var_set(1,5));
    E_XGY = X' * G_new * Y;

    [training_error1, norm_value1] = calculate_error(E_new, testing_Data, dump_matrix);
    [training_error2, norm_value2] = calculate_error(E_XGY, testing_Data, dump_matrix);
    training_error = min(training_error1, training_error2);
end