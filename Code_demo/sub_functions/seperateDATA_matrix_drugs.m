function [training_E, testing_E, dump_matrix] = seperateDATA_matrix_drugs(E, dump_percent, previous_dump_matrix)
[n1, n2] = size(E); % dimension of the tensor

drugs = floor(n2/11);

dump_matrix = false(n1, n2);
for i = 1:drugs
    temp = rand(n1, 1);
    temp = temp < dump_percent;
    temp = repmat(temp, 1, 11);
%     dump_position = permute(dump_position,[1 3 2]);
    dump_matrix(:, 11 * i-10:11 * i) = temp;
end
clear temp i drugs

training_E = E;
training_E(dump_matrix) = 0;
testing_E = zeros(n1, n2);
testing_E(dump_matrix) = E(dump_matrix);
% dump_matrix = dump_matrix(previous_dump_matrix);

end


