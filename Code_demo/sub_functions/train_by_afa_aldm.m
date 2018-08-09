function [result_table, best_E, best_G, X_shuffled, X_shuffled] = train_by_afa_aldm(method,X,Y,E_org,train_delete_rate, valid_delete_rate,lambda_seq_E,lambda_seq_G,seq_zeta)
    %% initialize parameters
    tic;
    maxIter = 3000;
    eps = 10e-6;
    
    %% shuffle matrix
    [X_shuffled,E_shuffled,shuffle_index] = shuffle_matrix(X,E_org);
    
    %% split matrix
    [train_matrix_org, valid_matrix_org, test_matrix_org] = split_3_matrix(E_shuffled); 
                
    %% generate train data
    [train_matrix, dumped_train_matrix] = generate_training_data(train_matrix_org, train_delete_rate);
    [valid_matrix, dumped_valid_matrix] = generate_training_data(train_matrix_org, valid_delete_rate);
    E = [train_matrix; valid_matrix; test_matrix_org];
    clear train_matrix valid_matrix
    toc;
    
    %% test with different parameters
    counter = 1;
    result_table = [];
    for i=1:length(lambda_seq_E)
        lambda_E =lambda_seq_E(i);
        for j=1:length(lambda_seq_G)
            lambda_G =lambda_seq_G(j);
            for k=1:length(seq_zeta)
                zeta = seq_zeta(k);
                
                %% run AFA_ALDM algorithm
                tic;
                if method == 'AFA_ALDM'                
                    [G_new, E_new] = AFA_ALDM(X, Y, E, lambda_E, lambda_G, maxIter, eps, zeta);
                    elapsedTime = toc;
                    E_XGY = X' * G_new * Y;
                    clear G_new
                    
                    [train_matrix_XGY, valid_matrix_XGY, ~] = split_3_matrix(E_XGY); 
                    [train_error_rate_XGY, train_norm_XGY] = calculate_error(train_matrix_org, train_matrix_XGY, dumped_train_matrix, train_delete_rate);
                    [valid_error_rate_XGY, valid_norm_XGY] = calculate_error(valid_matrix_org, valid_matrix_XGY, dumped_valid_matrix, valid_delete_rate);
                    clear train_matrix_XGY valid_matrix_XGY E_XGY

                    [train_matrix_new, valid_matrix_new, ~] = split_3_matrix(E_new); 
                    [train_error_rate_new, train_norm_new] = calculate_error(train_matrix_org, train_matrix_new, dumped_train_matrix, train_delete_rate);
                    [valid_error_rate_new, valid_norm_new] = calculate_error(valid_matrix_org, valid_matrix_new, dumped_valid_matrix, valid_delete_rate);
                    clear train_matrix_new valid_matrix_new E_new

                    result_table(counter,:) = [counter, lambda_E, lambda_G, zeta, elapsedTime, train_error_rate_XGY, train_error_rate_new, valid_error_rate_XGY, valid_error_rate_new, train_norm_XGY, train_norm_new, valid_norm_XGY, valid_norm_new];
                    clear elapsedTime train_error_rate_XGY train_error_rate_new valid_error_rate_XGY valid_error_rate_new train_norm_XGY train_norm_new valid_norm_XGY valid_norm_new
                    
                elseif method == 'dirtyIMC'
                    [UU SS VV U S V] = dirtyIMC(E, X', Y', lambda_E, lambda_G, 5, zeta); 
                    elapsedTime = toc;
                    G_new = UU*SS*VV';
                    E_new = X'*G_new*Y + U*S*V';
                    clear UU SS VV U S V G_new
                    
                    [train_matrix_new, valid_matrix_new, ~] = split_3_matrix(E_new); 
                    [train_error_rate_new, train_norm_new] = calculate_error(train_matrix_org, train_matrix_new, dumped_train_matrix, train_delete_rate);
                    [valid_error_rate_new, valid_norm_new] = calculate_error(valid_matrix_org, valid_matrix_new, dumped_valid_matrix, valid_delete_rate);
                    clear train_matrix_new valid_matrix_new E_new
                    
                    result_table(counter,:) = [counter, lambda_E, lambda_G, zeta, elapsedTime, train_error_rate_new, valid_error_rate_new, train_norm_new, valid_norm_new];
                    clear elapsedTime train_error_rate_new valid_error_rate_new train_norm_new valid_norm_new
                end
               
                counter = counter + 1;          
            end
        end
    end
end
