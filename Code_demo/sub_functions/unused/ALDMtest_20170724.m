%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
% The code is a demo of the AFA_ALDM method on movieLens.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('TCSIM_ALDM'))
addpath(genpath('code'));


%% read data. Notice we use 0 represent missing data in our matrix
% E = X'GY
% E is n by m matrix
% X is feature p by n matrix
% Y is feature q by m matrix

% E = X'GY
% E: 6678*22

% X:6678*389
% G:389*22
% Y:22*22



% Variable initialization
X = transpose(csvread('../data/processed data/sorted_combined_patient_feature_matrix_without_all_empty.csv'));
Y = csvread('../data/processed data/unit_matrix_with_dimension_22.csv');
E_org = csvread('../data/processed data/sorted_patient_symp_matr_with_o_all_emp_corrected.csv');

X = normc(X');
X = X';
X_100 = X(:, 1:100);

counter = 0;
E_org_full_observe = E_org;
i = 1;

while counter < 100
    if E_org_full_observe(i,1) == 0 | E_org_full_observe(i,12) == 0
        E_org_full_observe(i, :) = [];
    else
        counter = counter + 1;
        i = i + 1;
    end
end
clear counter;
clear i;

E_org_full_observe_100 = E_org_full_observe(1:100, :);
E_org_train_100 = E_org_full_observe_100;

random_index = sort(randperm(100, 50));
random_index_coc_or_opi = [];
for k=1:length(random_index)
    if rand >= 0.5
        for i=1:11
            E_org_train_100(random_index(k),i) = 0;
        end
        random_index_coc_or_opi(end+1) = 0;
    else
        for i=12:22
            E_org_train_100(random_index(k),i) = 0;
        end
        random_index_coc_or_opi(end+1) = 1;
    end
end

%% run AFA_ALDM model
maxIter = 3000;
eps = 10e-6;
zeta = 10e5;

result_norm = []
lambda_array = [0.01,0.1,1,10,100,1000]
counter = 1
for i=1:6
    for j=1:6
        lambda_E =lambda_array(i);
        lambda_G =lambda_array(j);

        tic;
        [G_new, E_new] = AFA_ALDM(X_100, Y, E_org_train_100,  lambda_E, lambda_G, maxIter, eps, zeta);
        elapsedTime = toc;

        norm_matrix = [];
        norm_matrix_rounded = [];
        for k=1:length(random_index)
            for l=1:11
                if random_index_coc_or_opi(k) == 0
                    norm_matrix(k,l) = E_org_full_observe_100(random_index(k), l) - E_new(random_index(k), l);
                    if E_new(random_index(k), l) > 0
                        norm_matrix_rounded(k,l) = E_org_full_observe_100(random_index(k), l) - 1;
                    else
                        norm_matrix_rounded(k,l) = E_org_full_observe_100(random_index(k), l) + 1;
                    end
                else
                    norm_matrix(k,l) = E_org_full_observe_100(random_index(k), l+11) - E_new(random_index(k), l+11);
                    if E_new(random_index(k), l+11) > 0
                        norm_matrix_rounded(k,l) = E_org_full_observe_100(random_index(k), l+11) - 1;
                    else
                        norm_matrix_rounded(k,l) = E_org_full_observe_100(random_index(k), l+11) + 1;
                    end
                end
            end
        end
        result_norm(counter,1) = counter;
        result_norm(counter,2) = lambda_E;
        result_norm(counter,3) = lambda_G;
        result_norm(counter,4) = elapsedTime;
        result_norm(counter,5) = norm(norm_matrix, 'fro') * 0.25 / (22*100);
        result_norm(counter,6) = norm(norm_matrix_rounded, 'fro') * 0.25 / (22*100);
        counter = counter + 1;
    end
end

csvwrite('test_result_for_100_patient.csv',result_norm)

clear lambda_E, lambda_G, i, j, k, l, lambda_array, random_index, random_index_coc_or_opi
