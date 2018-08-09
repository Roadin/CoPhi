%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
% The code is a demo of the Maxide[2] method on multi-label data Yahoo Arts
% [1] when the percentage of observed training instances' label assignments
% is 10%.
% 
% [1] N. Ueda, K. Saito. Parametric mixture models for multi-label text. 
% In: NIPS, 2002, 721-728.
% [2] M. Xu, R. Jin, Z-H. Zhou. Speedup matrix completion with side
% information: application to multi-label data. In: NIPS, 2013.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% STEP 0: Load Data, Set Default Value and Randomization Seed
addpath('Maxide');
dataname='Arts';
lambda=2^2;
percent=0.1;

data=['data/' dataname '.mat'];
load(data);

U=[train_data;test_data];
hatT=[train_target;test_target];
[n,q]=size(hatT);


trial=1;
s = RandStream.create('mt19937ar','seed',trial);
%For MATLAB version 2013a or later, using
RandStream.setGlobalStream(s);
%For MATLAB version lower than 2013a, using
%RandStream.setDefaultStream(s);

max_iter=100;
%% STEP 1: Generate Omega
train_num=round(n*0.9);

%Generate 90% training and 10% testing
obrT=zeros(size(hatT));
indexperm=randperm(n);
train_index=indexperm(1,1:train_num);
test_index=indexperm(1,train_num+1:n);
remainT=hatT(train_index,:);

%Generate the percent% randomly observed entries in the training data
for iii=1:q
    positive_index=find(remainT(:,iii)>0);
    positive_number=length(positive_index);
    positive_random=randperm(positive_number);
    positive_select=positive_index(positive_random(1,1:ceil(positive_number*percent)),1);
    negative_index=find(remainT(:,iii)<=0);
    negative_number=length(negative_index);
    negative_random=randperm(negative_number);
    negative_select=negative_index(negative_random(1,1:ceil(negative_number*percent)),1);
    obrT(train_index(1,positive_select),iii)=1;
    obrT(train_index(1,negative_select),iii)=1;
end

%% STEP 2: Recover M
disp('---------Matrix Recovery: ---------');
if  min(min(hatT))==-1
    hatT=(hatT+1)/2;
end
A=U;
M_Omega=hatT.*obrT;
Omega_linear=find(obrT);
B=eye(size(M_Omega,2));
[T,telapsed_side]=Maxide(M_Omega,Omega_linear,A,B,lambda,max_iter);

%% STEP 3: Evaluation
disp('---------Recovery Result: ---------');
results=PerformanceMeasure(T,hatT,test_index);
disp(['The average precision for the test instances is: ' num2str(results(1,1))]);
disp(['The average precision for the whole matrix is: ' num2str(results(1,2))]);
disp(['The training time is ' num2str(telapsed_side) 'seconds.'])