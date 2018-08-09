%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
% The code is a demo of the AFA_ALDM method on movieLens.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('TCSIM_ALDM'));
addpath(genpath('TCSIM_ALDM\code'));
addpath(genpath('compared methods\dirtyIMC_code 2'));
addpath(genpath('compared methods\imc-matlab_from_dhillion'));
addpath(genpath('compared methods\Maxide'));
addpath(genpath('sub functions'));


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

%% Variable initialization
X = csvread('../data/processed data/sorted_combined_patient_feature_matrix_without_all_empty.csv')';
Y = csvread('../data/processed data/unit_matrix_with_dimension_22.csv');
E_org = csvread('../data/processed data/sorted_patient_symp_matr_with_o_all_emp_corrected.csv');
X = normc(X');


%% generate full observe data
X_full_observe = zeros(size(X));
E_full_observe = zeros(size(E_org));
rows = size(E_org, 1);
counter = 1;
for i=1:rows
    if E_org(i,1) ~= 0 && E_org(i,12) ~= 0
        X_full_observe(counter, :) = X(i,:);
        E_full_observe(counter, :) = E_org(i,:);
        counter = counter +1;
    end
end
X_full_observe(counter:end, :) = [];
E_full_observe(counter:end, :) = [];

clear rows columns i counter
X_full_observe = X_full_observe';


%% create diagonal matrix Y_2
Y_2 = eye(22,22);
for i=1:11
    Y_2(i, i+11) = 1;
    Y_2(i+11, i) = 1;
end
clear i


%% generate first 30 sample data
X_full_observe_30 = X_full_observe(:, 1:30);
E_full_observe_30 = E_full_observe(1:30, :);


%% generate first 300 sample data
X_full_observe_400 = X_full_observe(:, 1:400);
E_full_observe_400 = E_full_observe(1:400, :);

%% set parameters for AFA_ALDM
% method = 'AFA_ALDM';
% lambda_E_set= [0.1,1,10,100,1000];
% lambda_G_set = [0.1,1,5,10];
% zeta_set = [1e2,500,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10];

% lambda_E_set = [0.1,1,10];
% lambda_G_set = [0.1];
% zeta_set = [1e9];

%% set parameters for dirtyIMC
% method = 'dirtyIMC';
% 
% lambda_set = [0.1,1,10,100,1000];
% lambda_1_set = [0.1,1,10,100,1000];
% iteration = [500];

%% set parameters for IMC
% method = 'IMC';
% k_set = [1,5,10];
% lambda_set = [0.1,1,10,100,1000];
% iteration = [500];

%% set parameters for Maxide
method = 'Maxide';
lambda_set = [0.1,1,10,100,1000];

%% run AFA_ALDM model
Maxiter = 5;
for k=5:5
    best_array = [];
    start_time = cputime;
    percentage_deleted = 0.1*k;                     
    for i=1:Maxiter
%         [result_table, best_norm, best_E, best???_G] = test_case('AFA_ALDM',X_full_observe,Y, E_full_observe, percentage_deleted,0.5, lambda_E_set, lambda_G_set, zeta_set);
        [result_table, best_norm, best_E, best_G] = test_case('dirtyIMC',X_full_observe,Y_2, E_full_observe, percentage_deleted,0.5, lambda_set, lambda_1_set, iteration);
%         [result_table, best_norm, best_E, best_G] = test_case('IMC',X_full_observe,Y_2, E_full_observe, percentage_deleted,0.5, k_set, lambda_set, iteration);
%         [result_table, best_norm, best_E, best_G] = test_case('Maxide',X_full_observe,Y, E_full_observe, percentage_deleted,0.5, lambda_set, [1], [1]);
        

        for j=1:size(result_table,1)
            if strcmp(method, 'AFA_ALDM')
                if min(result_table(j, 10:11)) == best_norm
                    best_array(end+1, :) = result_table(j,2:end);
                end
            else
                if result_table(j, 8) == best_norm
                    best_array(end+1, :) = result_table(j,2:end);
                end
            end
        end
        csvwrite(strcat(method,'_Y2_', num2str(percentage_deleted), '_', int2str(i), '.CSV'),result_table);
        csvwrite(strcat(method,'_Y2_', num2str(percentage_deleted), '_', int2str(i), '_E_.CSV'),best_E);
        csvwrite(strcat(method,'_Y2_', num2str(percentage_deleted), '_', int2str(i), '_G_.CSV'),best_G);
    end
    end_time = cputime - start_time;
    best_array(end+1,:) = mean(best_array);
    csvwrite(strcat(method,'_Y2_', num2str(percentage_deleted), '_summary_.CSV'),best_array);
end

%         figure
%         surf(best_E, 'edgecolor', 'none'); view(0,90); axis tight;
%         xlabel('x');
%         ylabel('y');

clear 