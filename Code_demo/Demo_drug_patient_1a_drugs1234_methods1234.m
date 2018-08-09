%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% The code is a demo of all matrix completion methods
%   Y_select        - side information y. Available options are:
%       'Ya'        - identity matrix
%       'Yb'        - diagonal matrix
%       'Yc'        - similarity matrix
%           'Ycc'   - cosine similarity
%           'Yci'   - infinite norm similarity
%           'Ycn'   - n norm similarity, usually we choose n=13
%       'Yd'        - association matrix
%   start_p, end_p  - dump_percent.
%                     Ranged from 0.1 to 0.9 with 0.1 interval.
%   test            - small sample test. Available optiosare 30 or 300                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load path and set result path
addpath(genpath('sub_functions'));
addpath(genpath('methods/LADMM_sto'));
addpath(genpath('methods/dirty_IMC'));
addpath(genpath('methods/IMC'));
addpath(genpath('methods/Maxide'));

data_path = strcat(pwd, '/data/processed_data/');
data_name = ['coc'; 'opi'; 'alc';'can'];
result_path = 'result_demo';
if ~(exist(result_path, 'file') == 7)
    mkdir(sprintf(result_path));
end
result_path = strcat(pwd, '/', result_path,'/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters configuration
dump_percent = [0.1, 0.3, 0.5, 0.7, 0.9];
repeat = 3;
Y_select = 'Ya';
% Select a small part of the data set to run a sample test.
% Disable it by sample_test = 0;
sample_test = 0;
% sample_test = 30;

%% Set parameters for each method
% Method 1: LADMM_sto
% Method 2: LADMM
% Method 3: Maxide
% Method 4: IMC
% Method 5: dirtyIMC
% Method = [1,2,3,4];
Method = [1];
Methods = size(Method,2);

% set parameters for LADMM_sto
LADMM_sto.lambda_E_set = [100,1000];
LADMM_sto.lambda_G_set = [0.1,1];
LADMM_sto.zeta_set = [10e7,10e8, 10e9];
LADMM_sto.mult_set = [8];
LADMM_sto.beta_0_set = [4*1e-1,1e-4,1e-6,1e-8];

% set parameters for LADMM
LADMM.lambda_E_set = [100,1000];
LADMM.lambda_G_set = [0.1,1];
LADMM.zeta_set = [10e7,10e8, 10e9];
LADMM.mult_set = [8];
LADMM.beta_0_set = [4*1e-1,1e-4,1e-6,1e-8];

% set parameters for dirtyIMC_params
dirtyIMC_params.lambda_set = [1];
dirtyIMC_params.lambda_1_set = [10];
dirtyIMC_params.iteration = [5];

% set parameters for IMC_params
IMC_params.k_set = [10,100,1000];
IMC_params.lambda_set = [10,100,1000];
IMC_params.iteration = [100];

% set parameters for Maxide
Maxide_params.lambda_set = [100,1000];

% % set parameters for inexact_alm_rpca
% IAR.lambda_set = [0,0.1,1,10,100,1000];
% IAR.eps_set = [1e-5, 1e-6,1e-7,1e-8,1e-9,1e-10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data from files

% Read drug_sympton data for each drug
drugs = 4;
F_full_observe = zeros(3131,11*drugs);
for i=1:drugs
        F_filename = strcat(data_path, 'F_', data_name(i,:), '-data');
        F_full_observe(:,11*i-10:11*i) = csvread(strcat(F_filename, '_full_observe_matched.csv'),0,1);
        
end
clear F_filename
F_full_observe(F_full_observe == 0) = -1;

% Read X_full_observe
X_full_observe = csvread(strcat(data_path,'X_patident_feature_full_observe_matched.csv'),0 , 1);
X_full_observe = normc(X_full_observe)';
%% Sample test (Select a small part of the data set to run a sample test)
if sample_test > 0
    %% generate first 400 sample data
    X_full_observe = X_full_observe(:, 1:sample_test);
    F_full_observe = F_full_observe(1:sample_test, :);
end

%% Y generation
Y = y_generation(Y_select, X_full_observe, F_full_observe);

%% run test case
input_parameters.result_path = result_path;
input_parameters.Method = Method;
input_parameters.Methods = Methods;
input_parameters.LADMM = LADMM;
input_parameters.LADMM_sto = LADMM_sto;
input_parameters.dirtyIMC_params = dirtyIMC_params;
input_parameters.IMC_params = IMC_params;
input_parameters.Maxide_params = Maxide_params;

input_parameters.F = F_full_observe;
input_parameters.X= X_full_observe;
input_parameters.Y = Y;
input_parameters.dump_percent = dump_percent;
input_parameters.repeat =repeat;

[Error_mean, Error_var, Elapsed_time] = relative_Error_var(input_parameters);        

%% Write result to files
csvwrite(strcat(result_path, 'Error_mean.csv'), Error_mean);
csvwrite(strcat(result_path, 'Error_var.csv'), Error_var);
csvwrite(strcat(result_path, 'Elapsed_time.csv'), Elapsed_time);
