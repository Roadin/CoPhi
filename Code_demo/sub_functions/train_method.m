function     [Elapsed_time, norm_E, norm_G, best_G] = train_method(method, X, Y, F, dump_position, E_org, G, variable_1_set,variable_2_set,variable_3_set)

    maxIter = 1000;
    eps = 10e-6;
    zeta = 10e6;
    L1 = length(variable_1_set);
    L2 = length(variable_2_set);
    L3 = length(variable_3_set);
    
    temp_table = zeros(L1*L2*L3, 3);
    G_cell = cell(L1*L2*L3,1);
    
    for i=1:L1
        variable_1 =variable_1_set(i);
        for j=1:L2
            variable_2 =variable_2_set(j);
            for k=1:L3
                variable_3 =variable_3_set(k);
                if strcmp(method, 'AFA_ALDM')           
                    [G_new_dump, E_new_dump,Elapsed_time] = AFA_ALDM(X', Y, F, variable_1, variable_2, maxIter, eps, zeta);
                    
                elseif strcmp(method, 'AFA_ALDM_sto')
                    [G_new_dump, E_new_dump, Elapsed_time] =  AFA_ALDM_sto(X', Y, F,  variable_1, variable_2, maxIter, eps,zeta,variable_3); 

                    
                elseif strcmp(method, 'dirtyIMC')
                    tic;
                    [UU, SS, VV, U, S, V] = dirtyIMC(F, X, Y', variable_1, variable_2, 100, 0, eps); 
                    Elapsed_time = toc;
                    G_new_dump = UU*SS*VV';
                    clear UU SS VV
                    E_new_dump = X*G_new_dump*Y + U*S*V';
                    clear U S V

                elseif strcmp(method, 'IMC')
                    [W, H, Elapsed_time] = IMC(F, X, Y', variable_1, variable_2, maxIter, [], []);
                    G_new_dump = W'*H;
                    clear W H
                    E_new_dump = X*G_new_dump*Y;

                elseif strcmp(method, 'Maxide')
                    F_vector = F(:);
                    F_linear = [];
                    for m=1:size(F_vector,1)
                        if F_vector(m,1)~=0
                            F_linear(end+1) = m;
                        end
                    end
                    clear F_vector
                    tic;
                    [E_new_dump,G_new_dump] = Maxide(F, F_linear, X, Y', variable_1, maxIter);
                    Elapsed_time = toc;
                    clear E_linear
                end
                current_line = (i-1)* L2*L3+ (j-1)* L3+ k;
                temp_table(current_line, 1) = Elapsed_time;
                temp_table(current_line, 2) = norm( E_new_dump( dump_position) - E_org( dump_position))/norm(E_org( dump_position));
                temp_table(current_line, 3) = norm(G_new_dump - G) / norm(G);
                G_cell{current_line} = G_new_dump;
            end
        end
    end
    best_norm_G = min(temp_table(:,3));
    for i=1:size(temp_table, 1)
        if temp_table(i,3) == best_norm_G
            Elapsed_time = temp_table(i,1);
            norm_E = temp_table(i,2);
            norm_G = temp_table(i,3);
            best_G = G_cell{current_line};
        end
    end
end