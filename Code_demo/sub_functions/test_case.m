function [Error_mean, Error_var, Elapsed_time] = test_case(params)
    
    params.Methods = Methods
    dump_percent = dumpt_percent.dump_percent;
    Error_mean = zeros(Methods, length(dump_percent));
    Error_var = Error_mean;
    Elapsed_time = Error_mean;
    result_path = params.result_path;
    
    F = input_parameters.F;
    X = input_parameters.X;
    Y = input_parameters.Y;
    
    
    for i = 1:length( dump_percent)
        fprintf('Outestloop: dump_percent = %d\n', dump_percent(i));
        relative_Error_mean = 0;
        relative_Error_var = 0;
        relative_Elapsed_time = 0;
        G_new_sto = 0;
        G_new = 0;
        
        Error_mean(:, i) = relative_Error_mean;
        Error_var(:, i) = relative_Error_var;
        Elapsed_time(:, i) = relative_Elapsed_time;
        
        csvwrite(strcat(result_path,'current_Error_mean.csv'), Error_mean);
        csvwrite(strcat(result_path,'current_Error_var.csv'), Error_var);
        csvwrite(strcat(result_path,'current_elapsed_table.csv'), elapsed_table);
        csvwrite(strcat(result_path,'C_new_',num2str(dump_percent(i)),'.csv'), C_new);
        csvwrite(strcat(result_path,'C_new_sto_',num2str(dump_percent(i)),'.csv'), C_new_sto);  
    end
    
    % initialize parameters
    maxIter = 3000;
    eps = 10e-6;
    
    method = parameters.method;
    X = parameters.X;
    Y = parameters.Y;
    E_org = parameters.F;
    training_rate = parameters.t_rate;
    validation_rate = parameters.v_rate;
    variable_matrix = parameters.variable_matrix;
    variable_size = parameters.variable_size;
    [m, n] = size( E_org );

    
    %% shuffle matrix
    [X_shuffled,E_shuffled,shuffle_index_1] = shuffle_matrix(X,E_org);
%   [X_back,E_back]=unshuffle_matrix(X_s2,E_s2,shuffle_index_2);
                
    %% generate train data
    delete_rate = training_rate+(1-training_rate) * validation_rate;
    remove_index = randperm(m,floor(delete_rate*m))';
    remove_index(:,2) = randi([0,1],size(remove_index,1),1);
    clear delete_rate
    
    train_index = sortrows(remove_index(1:floor(training_rate*m),:),1);
    valid_index = sortrows(remove_index(floor(training_rate*m)+1:end,:),1);
    
    train_full_index = false(m,n);
    valid_full_index = false(m,n);
    train_full_index(train_index(:,1),1:11) = repmat((train_index(:,2) == 0),1,11);
    train_full_index(train_index(:,1),12:22) = repmat((train_index(:,2) == 1),1,11);
    valid_full_index(valid_index(:,1),1:11) = repmat((valid_index(:,2) == 0),1,11);
    valid_full_index(valid_index(:,1),12:22) = repmat((valid_index(:,2) == 1),1,11);
    
    E_valid = E_shuffled;
    E_train(train_full_index) = 0;
    E_train = E_valid;
    E_train(valid_full_index) = 0;
    E = E_train;
    X_s2 = X_shuffled;
    %% test with different parameters
    counter = 1;
    best_norm = 0;
    best_E = [];
    best_G = [];
    
    %% accelerate calculation by pre_process some data    
    if strcmp(method, 'AFA_ALDM') || strcmp(method, 'AFA_ALDM_sto') 
        dump_percent = 0.1;

        %% p is features of matrix X; q is features of matrix Y
        p = size(X, 1);
        q = size(Y, 1);
        G = zeros(p, q);

        %% First, we use SVT  to find initial matrix
        temp = reshape( E, m * n, 1 );
        [Tau, col] = find( ~(temp == 0));
        data = temp( Tau );
        tau = 5*sqrt(m * n); 
        delta = 1;
        maxit = floor( 500/dump_percent);
        tol = 1e-1;
        [U, S, V, numiter] = SVT(size(E), Tau, data, tau, delta, maxit, tol);
        E_2 =  U*S*V';
        E_2(E~=0) = E(E~=0);
        A = kron(Y', X');
    end
    
    
    
    
    
    result_table = zeros(variable_size,17);
    iter=0;
    vm = variable_matrix;
    
    
    valid = 0;
    i = 1;
    while i <= variable_size                
        %% run AFA_ALDM algorithm and calculate error
        display(counter);
        if valid
            E = E_valid;
        else
            E = E_train;
        end
        
        start_time = cputime;
        
        
        %% run method AFA_ALDM
        if strcmp(method, 'AFA_ALDM')
            [G_new, E_new,~] = AFA_ALDM(X_s2, Y, E, vm(i,1), vm(i,2), maxIter, eps, 100);

            elapsedTime = cputime - start_time;
            clear start_time
            E_XGY = X_s2' * G_new * Y;
        
        %% run method AFA_ALDM_sto
        elseif strcmp(method, 'AFA_ALDM_sto')
%           [G_new, E_new,~] = AFA_ALDM_sto_modified(X_s2, Y, E, vm(i,1), vm(i,2), maxIter, eps, vm(i,3),vm(i,4),E_2,A,vm(i,5));
            [G_new, E_new,~,iter] = AFA_ALDM_sto(X_s2, Y, E, vm(i,1), vm(i,2), maxIter, eps, vm(i,3),vm(i,4));
            elapsedTime = cputime - start_time;
            clear start_time
            E_XGY = X_s2' * G_new * Y;

        %% run method inexact_alm_rpca
        elseif strcmp(method, 'inexact_alm_rpca')
            [E_XGY, E_new, ~] = inexact_alm_rpca(E, vm(i,1), vm(i,2), maxIter);
            elapsedTime = cputime - start_time;
            clear start_time
            G_new = E_XGY;

        %% run method dirtyIMC
        elseif strcmp(method, 'dirtyIMC')
            [UU, SS, VV, U, S, V] = dirtyIMC(E, X_s2', Y', vm(i,1), vm(i,2), vm(i,3), 0, eps); 
            elapsedTime = cputime - start_time;
            G_new = UU*SS*VV';
            clear UU SS VV
            E_new = X_s2'*G_new*Y + U*S*V';
            clear U S V
            E_XGY = E_new;
        
        %% run method IMC
        elseif strcmp(method, 'IMC')
            [W, H, ~] = IMC(E, X_s2', Y, vm(i,1), vm(i,2), vm(i,3), [],[]);
            elapsedTime = cputime - start_time;
            G_new = W'*H;
            clear W H
            E_new = X_s2'*G_new*Y;
            E_XGY = E_new;
        
        %% run method Maxide
        elseif strcmp(method, 'Maxide')
            E_vector = E(:);
            E_linear = [];
            for m=1:size(E_vector,1)
                if E_vector(m,1)~=0
                    E_linear(end+1) = m;
                end
            end
            clear E_vector
            [E_new,G_new] = Maxide(E, E_linear, X_s2', Y, vm(i,1),maxIter);
            elapsedTime = cputime - start_time;
            clear E_linear
            E_XGY = E_new;
            
        %% run method naive0
        elseif strcmp(method, 'naive')
            E_new = E;
            for n=1:size(E,1)
                if E(n,1) == 0
                    E_new(n,1:11) = E_new(n,12:22);
                elseif E(n,12) == 0
                    E_new(n,2:12) = E_new(n,1:11);
                end
            end
            elapsedTime = 0;
            G_new = [];
            E_XGY = E_new;
            
        %% run method naive1
        elseif strcmp(method, 'naive1')
            E_new = E;
            for n=1:size(E,1)
                if E(n,1) == 0
                    E_new(n,1:11) = 1;
                elseif E(n,12) == 0
                    E_new(n,2:12) = 1;
                end
            end
            elapsedTime = 0;
            G_new = [];
            E_XGY = E_new;
        end

        [train_error_rate_XGY, train_norm_XGY] = calculate_error(abs(E_XGY-E_shuffled), valid_full_index);
        [train_error_rate_new, train_norm_new] = calculate_error(abs(E_new-E_org), valid_full_index);

        current.norm = min(train_norm_XGY, train_norm_new);
        current.norm_XGY = train_norm_XGY;
        current.E_XGY = E_XGY;
        current.G_XGY = G_XGY;

        current.norm_new = train_norm_new;
        current.E_new = E_new;
        current.G_XGY = G_new;

        
        [best_norm, best_E, best_G] = find_best(best_norm, best_E, best_G,current_norm, current_E, current_G);               
        clear current_norm current_E current_G

        %% save result in table
        result_table(counter,:) = [iter, vm(i,1), vm(i,2), vm(i,3),vm(i,4),vm(i,5), elapsedTime,-1, train_error_rate_XGY, train_error_rate_new, test_error_rate_XGY, test_error_rate_new,-1, train_norm_XGY, train_norm_new, test_norm_XGY, test_norm_new];
        clear elapsedTime train_error_rate_XGY train_error_rate_new test_error_rate_XGY test_error_rate_new train_norm_XGY train_norm_new test_norm_XGY test_norm_new
        csvwrite('current.csv', result_table);
        counter = counter + 1;
        
        if valid
            i = i + 1;
            valid = 0
        else
            valid = 1;
        end
    end
    [~,best_E] = resumE_s2_matrix(X_s2,best_E,shuffle_index_2);
end
