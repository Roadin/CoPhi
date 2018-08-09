function [variable_matrix, variable_size] = generate_variable_permutation(variable_set_1,variable_set_2,variable_set_3,variable_set_4,variable_set_5)
    variable_size = length(variable_set_1) * length(variable_set_2) * length(variable_set_3) * length(variable_set_4) * length(variable_set_5);

    vs1 = variable_set_1;
    vs2 = variable_set_2;
    vs3 = variable_set_3;
    vs4 = variable_set_4;
    vs5 = variable_set_5;
    counter = 1;
    variable_matrix = zeros(variable_size,5);
    for i1 = 1:length(vs1)
        for i2 = 1:length(vs2)
            for i3 = 1:length(vs3)
                for i4 = 1:length(vs4)
                    for i5 = 1:length(vs5)
                        variable_matrix(counter,:) = [vs1(i1), vs2(i2), vs3(i3), vs4(i4), vs5(i5)];
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end