function [ training_error, lambda1, lambda2, k1, k2] = CV_dirtyimc( X, Y, trainingData, testing_Data, lambda_seq)


    %%10folder CV --- tune parameter for DirtyIMC algorithm
    %parameters

    training_error_matrix = zeros( length( lambda_seq), length( lambda_seq));


    for i = 1:length(lambda_seq)
        for j = 1:length( lambda_seq)
            lambda1 = lambda_seq( i );
            lambda2 = lambda_seq( j );

            relative_error = zeros( 1, 5);
            for m = 1 : 5
                 fprintf('         i=%d, j=%d, m=%d\n', i, j, m);
                 [data_after_dump, dump_position] = selectdata( trainingData, 0.1);
%                  [W, H, U, V] = Dirty_IMC(data_after_dump, X', Y', lambda1, lambda2, maxiter);
                 [UU, SS, VV, U, S, V] = dirtyIMC(data_after_dump, X', Y', lambda1, lambda2);
                 M = UU*SS*VV'; 
                 N = U*S*V'; 

                 E_new_dump = X'*M*Y + N;
                 relative_error( m ) =  norm( E_new_dump( dump_position) - trainingData( dump_position))/norm(trainingData( dump_position));
            end    
            training_error_matrix( i, j) = mean( relative_error);

        end
    end    
    index = find( training_error_matrix == min( min( training_error_matrix )));
    index = index(1);
    [I_row, I_col] = ind2sub(size( training_error_matrix), index);

    
    lambda1 = lambda_seq( I_row );
    lambda2 = lambda_seq( I_col );
    
    %train dirty_IMC model
    %[W, H, U, V] = Dirty_IMC(trainingData, X', Y', k1, k2, lambda1, lambda2, maxiter);
    [UU, SS, VV, U, S, V] = dirtyIMC(trainingData, X', Y', lambda1, lambda2);
     M = UU*SS*VV'; 
     N = U*S*V'; 

     E_new= X'*M*Y + N;
    Omega  = testing_Data>0;
    training_error =  norm( E_new( Omega) - testing_Data( Omega))/norm(testing_Data( Omega));
    k1 = rank( M );
    k2 = rank( N );


end
