clear all

% addpath(genpath('/home/users/taghia/BSDS/'))

% generate a random data 
nsubjs = 3;
dim = 4;
nsamps = 60;
for subj =1 : nsubjs
    data{subj} = randn(dim, nsamps);
end

%% train BSDS model using defualt setting 

max_states = 20; % or any other sensible values
group_model = BayesianSwitchingDynamicalSystems(data, max_states);

%% train BSDS model using default setting for large number of ROIs
%NOTE: In applications where  n_rois (dim) is large, you would need to set max_ldim to smaller values (max_ldim < dim) and effectively limit the model complexity. This can be done by setting max_ldim to account for certain variaitons in the data (e.g., 95%). 

% [max_ldim, ~] = decideOnLocalComplexityBasedOnEnergy(data, 0.95); 
% max_states = 15; 
% group_model = BayesianSwitchingDynamicalSystems(data, max_states, max_ldim);


%% train BSDS model using advanced setting

% max_ldim = size(data{1}, 1)-1;
% max_states = 10; 
% opt.n_iter = 100;
% opt.n_init_iter =2; 
% opt.tol = 1e-10;
% opt.noise = 0;
% opt.n_init_learning = 2;
% group_model = BayesianSwitchingDynamicalSystems(data, max_states, max_ldim, opt);

%% compute subject-level covariance matrices from the trained group model

% subj_model = compute_subject_level_stats(data, group_model, opt.n_iter);

%% compute autoregressive coefficients from the group BSDS model

% [stats_group, stats_subj] = vector_autoregressive_model(data, group_model, subj_model);


