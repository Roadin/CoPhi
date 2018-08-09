function [relative_Error_mean, relative_Error_variance, relative_Elapsed_time, G_new_sto, G_new] = relative_error_fixed_dump(params, dump_percent)
    Methods = params.Methods;
    repeat = params.repeat;
    relative_Error_matrix = zeros(repeat, Methods);
    relative_Elapsed_matrix = relative_Error_matrix;
    result_path = params.result_path;
    
    F = params.F;
    
    for i = 1: repeat
        fprintf( 'Second-layer-loop: repeated sampling = %d\n', i);
        work_dir = strcat(result_path, num2str(dump_percent), '/seperateDATA/DATA', num2str(i));
        mkdir(work_dir);

        [training_Data, testing_Data, dump_matrix] = seperateDATA_matrix_drugs(F, dump_percent);

        csvwrite( strcat(work_dir, '/training_Data.csv'), training_Data);
        csvwrite( strcat(work_dir, '/testing_Data.csv'), testing_Data);
        [result_table ,G_new,G_new_sto] = CV_for_ALL_fixed_dump_matrix(params, training_Data, testing_Data, dump_matrix, work_dir, dump_percent);
        relative_Error_matrix( i, :) = result_table(1,:);
        relative_Elapsed_matrix(i,:) = result_table(2,:);
    end

    relative_Error_mean = zeros( 1, Methods);
    relative_Error_variance = zeros( 1, Methods);
    elapsed = zeros(1,Methods);
    for i = 1:Methods
        relative_Error_mean( i ) = mean( relative_Error_matrix( :, i));
        relative_Error_variance( i ) = var(  relative_Error_matrix( :, i));
        relative_Elapsed_time(i) = mean(relative_Elapsed_matrix(:,i));
    end
end