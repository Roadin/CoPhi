function [result_table ,G_new,G_new_sto] = CV_for_ALL_fixed_dump(params, training_Data, testing_Data, dump_matrix, work_dir, dump_percent)
    
    % Error matrix, which is used to store relative error
    %first row represent AFA_aldm, Maxide, IMC_params and Dirty_IMC_params, 
    X = params.X;
    Y = params.Y;
    [m, n] = size(training_Data);
    G_new=0;
    G_new_sto = 0;
    Method = params.Method;
    Methods = params.Methods;
    result_table = zeros(2, Methods); % E_new performance
    [training_Data_inside, validating_Data, dump_matrix_validating] = seperateDATA_matrix_drugs(training_Data, 0.5);
%     X = X(1:30,:);
    dump_position = dump_matrix_validating;
    dump_position(dump_matrix) = 0;
   
    for k=1:Methods
        if Method(1,k) == 1
            % calculate relative error of AFA_aldm_sto
            fprintf('      train model for AFA_ALDM_sto\n');
            [training_error, Elapsed_time,G_new_sto] = CV_for_LADMM_sto(params.LADMM_sto, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix);

        elseif Method(1,k) == 2
            % calculate relative error of AFA_aldm_sto
            fprintf('      train model for AFA_ALDM_sto\n');
            [training_error, Elapsed_time,G_new] = CV_for_LADMM(params.LADMM, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix);

        elseif Method(1,k) == 3    
            % calculate relative error of Maxid
            fprintf('      train model for Maxid\n');
            
            [training_error, Elapsed_time] = CV_for_Maxide(params.Maxide_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix);
                
        elseif Method(1,k) == 4    
            %calculate relative error of IMC_params
            training_error = -1;
            fprintf('      train model for IMC_params\n');
            [training_error, Elapsed_time] = CV_for_IMC(params.IMC_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix);

        elseif Method(1,k) == 5
            fprintf('      train model for Dirty_IMC_params\n');
            training_error = -1;
            [training_error, Elapsed_time] = CV_for_dirtyIMC(params.dirtyIMC_params, X, Y, training_Data_inside, dump_position, training_Data, testing_Data, dump_matrix);

%     Detail_Parameters = [sprintf( 'AFA_ALDM - lambda_E: %d, lambda_G: %d\n',lambda_E, lambda_G),...
%         sprintf( 'Maxide - lambda: %d\n', lambda_maxide),...
%         sprintf( 'IMC_params - k: %d, lambda: %d\n', k_0, lambda_IMC_params),...
%         sprintf( 'DirtyIMC_params - k1: %d, k2: %d,  lambda1: %d, lambda2: %d\n', k1, k2, lambda1_dirty, lambda2_dirty )];
% 
% 
%     dlmwrite( sprintf('%sParametersDetail', work_dir), Detail_Parameters, 'delimiter', '')

%     dlmwrite( sprintf('%sresult_table', work_dir), result_table, 'delimiter', ',')
        end
        
        result_table(1,k) = training_error;
        result_table(2,k) = Elapsed_time;        

        csvwrite(strcat(work_dir, 'current_result.csv'), result_table);
    end
end