
addpath(genpath('sub functions'));
addpath(genpath('compared methods\ALM\PROPATH'));
E_coc_opi_full_observe = csvread('E_coc_opi_full_observe.csv');

E_alc_can_full_observe = csvread('alc_can_full_observe_new.csv');
E_alc_can_full_observe = E_alc_can_full_observe(:, 2:end);

% E_full_observe = E_coc_opi_full_observe;
% dumped_matrix = zeros(3441,22);

E_full_observe = E_alc_can_full_observe;
dumped_matrix = zeros(4500, 22);

E = E_full_observe; 
% E(1:500,1:11) = 0;
% dumped_matrix(1:500, 1:11) = 1;
% 
E(1:500,12:22) = 0;
dumped_matrix(1:500, 12:22) = 1;


%% initialize parameters
maxIter = 3000;
eps = 10e-6;

%% test with different parameters

lambda_set = [0,0.1,1,10,100,1000];
eps_set = [1e-5, 1e-6,1e-7,1e-8,1e-9,1e-10];



repeat_best = zeros(5,6);
for repeat_count = 1:5
    result_table = zeros(size(lambda_set, 2) * size(eps_set, 2), 6);
    counter = 1;
    for i=1:length(lambda_set)
        lambda =lambda_set(i);
        for j=1:length(eps_set)
            eps =eps_set(j);

            %% run AFA_ALDM algorithm and calculate error
            display(counter);
            start_time = cputime;
            [E_XGY, E_new, ~] = inexact_alm_rpca(E, lambda, eps, maxIter);
            elapsedTime = cputime - start_time;
            clear start_time
            [error_rate_XGY, ~] = calculate_error(E_full_observe, E_XGY, dumped_matrix, 250/size(E_full_observe, 1));
            [error_rate_new, ~] = calculate_error(E_full_observe, E_new, dumped_matrix, 250/size(E_full_observe, 1));

            %% save result in table
            result_table(counter,:) = [counter, lambda, eps, elapsedTime, error_rate_XGY, error_rate_new];

            counter = counter + 1;
        end
    end
    
    best_error_rate = min(result_table(:, 5));
    for i=1:size(lambda_set, 2) * size(eps_set, 2)
        if result_table(i,5) == best_error_rate
            repeat_best(repeat_count, :) = result_table(i,:);
        end
    end
end
repeat_best(end+1, :) = mean(repeat_best);


a = 1;
