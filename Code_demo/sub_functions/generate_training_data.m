function [training_data, dumped_matrix] = generate_training_data(E, percentage_deleted)
    
    [rows,columns] = size(E);
    deleted_rows = floor(rows*percentage_deleted * 2);
    dumped_matrix = zeros(rows, columns);
    training_data = E;
    random_index = sort(randperm(rows, deleted_rows));
    random_left_right = randi([0,1], 1, deleted_rows);

    for k=1:length(random_index)
        if random_left_right(k) == 0
            training_data(random_index(k),1:floor(columns/2)) = 0;
            dumped_matrix(random_index(k),1:floor(columns/2)) = 1;
        else
            training_data(random_index(k),floor(columns/2)+1:columns) = 0;
            dumped_matrix(random_index(k),floor(columns/2)+1:columns) = 1;
        end
    end