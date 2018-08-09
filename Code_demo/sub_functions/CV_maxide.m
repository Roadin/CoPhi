function [ training_error, lambda] = CV_maxide( X, Y, trainingData, testing_Data, lambda_seq)


    %%10folder CV --- tune parameter for maxide algorithm
    %parameters
    maxiter = 300;
    
    training_error_matrix = zeros( 1, length( lambda_seq));

    [n1, n2] = size( trainingData);
    for i = 1:length(lambda_seq)
        lambda= lambda_seq( i );

        relative_error = zeros( 1, 5);
        for j = 1 : 5
            fprintf('         i=%d, j=%d', i, j);
             temp = reshape( trainingData, n1*n2, 1);
             [data_after_dump, dump_position] = selectdata( trainingData, 0.1);
             temp( dump_position ) = 0;
             [Omega_dump, ~] = find( temp~=0);

             [E_new_dump,~] = Maxide(data_after_dump,Omega_dump,X, Y, lambda, maxiter);
             relative_error( j ) =  norm( E_new_dump( dump_position) - trainingData( dump_position))/norm(trainingData( dump_position));

        end    
        training_error_matrix( i) = mean( relative_error);

    end    
    index = find( training_error_matrix ==  min( training_error_matrix ));
    index = index(1);

    lambda = lambda_seq( index );
    
    %train maxide model
     temp = reshape( trainingData, n1*n2, 1);
     [Omega, ~] = find( temp~=0);

     [E_new,~] = Maxide(trainingData,Omega,X, Y, lambda, maxiter);
    
    Omega  = testing_Data>0;
    training_error =  norm( E_new( Omega) - testing_Data( Omega))/norm(testing_Data( Omega));


end
