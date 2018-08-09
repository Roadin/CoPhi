function [ training_error, lambda] = CV_imc( X, Y, trainingData, testing_Data, lambda_seq, k)


    %%10folder CV --- tune parameter for maxide algorithm
    %parameters
    maxiter = 100;
    
    training_error_matrix = zeros( 1, length( lambda_seq));

    d1 = size( X,1);
    d2 = size( Y,1);
    W0 = randn(d1, k);
    H0 = randn(d2, k);
    
    for i = 1:length(lambda_seq)
        lambda= lambda_seq( i );

        relative_error = zeros( 1, 5);
        for j = 1 : 5
            fprintf('         i=%d, j=%d\n', i, j);
             [data_after_dump, dump_position] = selectdata( trainingData, 0.1);
             [W, H] = IMC(data_after_dump, X', Y', k, lambda, maxiter, W0, H0);
             E_new_dump = X'*W'*H*Y;
             relative_error( j ) =  norm( E_new_dump( dump_position) - trainingData( dump_position))/norm(trainingData( dump_position));

        end    
        training_error_matrix( i) = mean( relative_error);

    end    
    index = find( training_error_matrix ==  min( training_error_matrix ));

    lambda = lambda_seq( index );
    lambda = lambda(1);
    
    %train maxide model
     [W, H] = IMC(trainingData, X', Y', k, lambda, maxiter, W0, H0);
     E_new = X'*W'*H*Y;
    
    Omega = testing_Data>0;
    training_error =  norm( E_new( Omega) - testing_Data( Omega))/norm(testing_Data( Omega));


end
