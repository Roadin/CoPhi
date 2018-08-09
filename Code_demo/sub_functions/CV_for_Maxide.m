function [training_error, Elapsed_time] = CV_for_Maxide(Maxide_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix)
   
    [var_set, var_size ] = generate_variable_permutation(Maxide_params.lambda_set, [-1],[-1],[-1],[-1]);
%     [X_U, ~, ~] = svd( X');
%     [ r_a, r_a2] = size( X );
%     if ~(r_a > r_a2)
%         X_U = X_U( :, 1:r_a);
%     end
%     [~, ~, Y_V] = svd( Y);
%     [ r_b, ~] = size( Y );
%     Y_V = Y_V( 1:r_b, : );
%     Y_V = Y_V';
    E_linear = find(training_Data_inside ~= 0);
    maxIter = 300;
%     relative_error_matrix_for_each_method = ones(var_size,1);
%     temp = training_Data_inside(:);
%     [Omega_dump, ~] = find( temp~=0);    
    
    
    for i=1:var_size
%         [E_new,~] = Maxide(training_Data_inside,Omega_dump,X_U, Y_V, var_set(i,1), maxIter);
        [E_new,~] = Maxide(training_Data_inside,E_linear,X', Y, var_set(i,1), maxIter);
        [relative_error_matrix_for_each_method(i,1),~] = calculate_error(E_new, training_Data, dump_position);    
    end
    best_relative_error = min(relative_error_matrix_for_each_method);
    temp = training_Data(:);
    [Omega_dump, ~] = find( temp~=0);
    validating_var_set = var_set(find(relative_error_matrix_for_each_method == best_relative_error),:);
    tic;
    E_linear = find(training_Data ~= 0);
    [E_new,~] = Maxide(training_Data,E_linear,X', Y, validating_var_set(1,1), maxIter);
    Elapsed_time = toc;
    [training_error, ~] = calculate_error(E_new, testing_Data, dump_matrix);    
end