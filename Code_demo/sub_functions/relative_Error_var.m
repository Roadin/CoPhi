function [Error_mean, Error_var, Elapsed_time] = relative_Error_var(params)
    
    Methods = params.Methods;
    dump_percent = params.dump_percent;
    Error_mean = zeros(Methods, length(dump_percent));
    Error_var = Error_mean;
    Elapsed_time = Error_mean;
    result_path = params.result_path;
    

    
    for i = 1:length( dump_percent)
        fprintf('Outestloop: dump_percent = %d\n', dump_percent(i));
        [relative_Error_mean, relative_Error_variance, relative_Elapsed_time, G_new_sto, G_new] = relative_error_fixed_dump(params, dump_percent(i));

        Error_mean(:, i) = relative_Error_mean;
        Error_var(:, i) = relative_Error_variance;
        Elapsed_time(:, i) = relative_Elapsed_time;
        
        csvwrite(strcat(result_path,'current_Error_mean.csv'), Error_mean);
        csvwrite(strcat(result_path,'current_Error_var.csv'), Error_var);
        csvwrite(strcat(result_path,'current_elapsed_table.csv'), Elapsed_time);
        csvwrite(strcat(result_path,'G_new_',num2str(dump_percent(i)),'.csv'), G_new);
        csvwrite(strcat(result_path,'G_new_sto_',num2str(dump_percent(i)),'.csv'), G_new_sto);  
    end
end
