function [result_table, best_norm, best_E, best_G] = test_case_sto(method,X,Y,E_org,train_delete_rate, test_delete_rate,variable_set_1,variable_set_2,variable_set_3)
    %% initialize parameters
    maxIter = 3000;
    eps = 10e-6;
    
    %% shuffle matrix
    [X_shuffled,E_shuffled,shuffle_index] = shuffle_matrix(X,E_org);
    
    %% split matrix
    [train_matrix_org, test_matrix_org, valid_matrix_org] = split_3_matrix(E_shuffled); 
                
    %% generate train data
    [train_matrix, dumped_train_matrix] = generate_training_data(train_matrix_org, train_delete_rate);
    [test_matrix, dumped_test_matrix] = generate_training_data(train_matrix_org, test_delete_rate);
    E = [train_matrix; test_matrix; valid_matrix_org];
    clear train_matrix test_matrix
    
    %% test with different parameters
    counter = 1;
    best_norm = 0;
    best_E = [];
    best_G = [];
    mult = 2;
    if strcmp(method, 'AFA_ALDM')
        result_table = zeros(size(variable_set_1,2)*size(variable_set_2,2)*size(variable_set_3,2),13);
    else
        result_table = zeros(size(variable_set_1,2)*size(variable_set_2,2)*size(variable_set_3,2),9);
    end
    
    for i=1:length(variable_set_1)
        variable_1 =variable_set_1(i);
        for j=1:length(variable_set_2)
            variable_2 =variable_set_2(j);
            for k=1:length(variable_set_3)
                variable_3 = variable_set_3(k);
                
                %% run AFA_ALDM algorithm and calculate error
                counter
                start_time = cputime;
                if strcmp(method, 'AFA_ALDM')           
                    [G_new, E_new] = AFA_ALDM_sto(X_shuffled, Y, E, variable_1, variable_2, maxIter, eps, variable_3, mult);
                    elapsedTime = cputime - start_time;
                    clear start_time
                    E_XGY = X_shuffled' * G_new * Y;
                    
                    [train_matrix_XGY, test_matrix_XGY, ~] = split_3_matrix(E_XGY); 
                    [train_error_rate_XGY, train_norm_XGY] = calculate_error(train_matrix_org, train_matrix_XGY, dumped_train_matrix, train_delete_rate);
                    [test_error_rate_XGY, test_norm_XGY] = calculate_error(test_matrix_org, test_matrix_XGY, dumped_test_matrix, test_delete_rate);
                    clear train_matrix_XGY test_matrix_XGY

                    [train_matrix_new, test_matrix_new, ~] = split_3_matrix(E_new); 
                    [train_error_rate_new, train_norm_new] = calculate_error(train_matrix_org, train_matrix_new, dumped_train_matrix, train_delete_rate);
                    [test_error_rate_new, test_norm_new] = calculate_error(test_matrix_org, test_matrix_new, dumped_test_matrix, test_delete_rate);
                    clear train_matrix_new test_matrix_new
                    
                    current_norm = [train_norm_XGY, train_norm_new];
                    current_E = cell(1,2);
                    current_G = cell(1,2);
                    current_E{1} = E_XGY;
                    current_E{2} = E_new;
                    current_G{1} = G_new;
                    current_G{2} = [];
                else
                    if strcmp(method, 'dirtyIMC')
                        [UU, SS, VV, U, S, V] = dirtyIMC(E, X_shuffled', Y', variable_1, variable_2, variable_3, 0, eps); 
                        elapsedTime = cputime - start_time;
                        G_new = UU*SS*VV';
                        clear UU SS VV
                        E_new = X_shuffled'*G_new*Y + U*S*V';
                        clear U S V

                    elseif strcmp(method, 'IMC')
                        [W, H, ~] = IMC(E, X_shuffled', Y, variable_1, variable_2, variable_3, [], []);
                        elapsedTime = cputime - start_time;
                        G_new = W'*H;
                        clear W H
                        E_new = X_shuffled'*G_new*Y;

                    elseif strcmp(method, 'Maxide')
                        E_vector = E(:);
                        E_linear = [];
                        for m=1:size(E_vector,1)
                            if E_vector(m,1)~=0
                                E_linear(end+1) = m;
                            end
                        end
                        clear E_vector
                        [E_new,G_new] = Maxide(E, E_linear, X_shuffled', Y, variable_1,maxIter);
                        elapsedTime = cputime - start_time;
                        clear E_linear
                    end

                    [train_matrix_new, test_matrix_new, ~] = split_3_matrix(E_new); 
                    [train_error_rate_new, train_norm_new] = calculate_error(train_matrix_org, train_matrix_new, dumped_train_matrix, train_delete_rate);
                    [test_error_rate_new, test_norm_new] = calculate_error(test_matrix_org, test_matrix_new, dumped_test_matrix, test_delete_rate);
                    clear train_matrix_new test_matrix_new

                    current_norm = train_norm_new;
                    current_E = cell(1,1);
                    current_G = cell(1,1);
                    current_E{1} = E_new;
                    current_G{1} = G_new;
%                 else
%                     err('no method found :(');
                end
                
                [best_norm, best_E, best_G] = find_best(best_norm, best_E, best_G,current_norm, current_E, current_G);               
                clear current_norm current_E current_G
                
                %% save result in table
                if strcmp(method, 'AFA_ALDM')
                    result_table(counter,:) = [counter, variable_1, variable_2, variable_3, elapsedTime, train_error_rate_XGY, train_error_rate_new, test_error_rate_XGY, test_error_rate_new, train_norm_XGY, train_norm_new, test_norm_XGY, test_norm_new];
                    clear elapsedTime train_error_rate_XGY train_error_rate_new test_error_rate_XGY test_error_rate_new train_norm_XGY train_norm_new test_norm_XGY test_norm_new

                else
                    result_table(counter,:) = [counter, variable_1, variable_2, variable_3, elapsedTime, train_error_rate_new, test_error_rate_new, train_norm_new, test_norm_new];
                    clear elapsedTime train_error_rate_new test_error_rate_new train_norm_new test_norm_new
                end
                
                counter = counter + 1;
                pause(5);
            end
        end
    end
    [~,best_E] = resume_shuffled_matrix(X_shuffled,best_E,shuffle_index);
end
