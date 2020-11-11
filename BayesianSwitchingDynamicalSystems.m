

function model = BayesianSwitchingDynamicalSystems(data, max_nstates, max_ldim, opt)

%Inputs:
%	data                             : {data_subj1, data_subj2,...,data_subjS}, data_subj1 is a D-by-N matrix where D is dimension  and N is the number of samples
%	max_nstates                 : maximum number of states (set to a value greater than the expected number of states)
%	max_ldim                      : bound on the latent space dimensionality. default: max_ldim = D-1. max_ldim must be smaller  than dimensionality of data
%	opt.n_iter                      : number of iteration of the variational inference
%   opt.i_init_iter                 : number of iterations for the initial learning
%   opt.tol                          : tolerance (1e-3)
%	opt.noise                      : type of measurement noise, can be either 1 or 0
%							                 - if opt.noise=1, noise variance the same across dimensions
%							                 - if opt.noise=0 (default), noise variance is allowed to vary across dimensions
%	opt.n_init_learning :    train a number of models using random initialization and choosing the best initial model (default: 10)

%Outputs:
%	  model.net 					                                : posterior parameters over all model parameters  (advanced usage)
%	  model.estimated_mean 		                            : estimated mean of each state (group-level)
%	  model.estimated_covariance                          : estimated covariance of each state (group-level)
%	  model.posterior_probabilities		                    : posterior probabilities of states across time for each subject q(Z)
%	  model.joint_posterior_probabilities	              : joint posterior probabilities of states across time for each subject q(Z_n-1, Z)
%	  model.state_transition_probabilities	             : HMM estimated transition probabilities across states
%     model.temporal_evolution_of_states               : estimated temporal evolution of states
%     model.fractional_occupancy_group_wise         : occupancy rate of each state computed in group sense
%     model.mean_lifetime_group_wise                   : mean life time of each state computed in group sense
%     model.id_of_dominant_states_group_wise       : id of dominant states group wise
%     model.fractional_occupancy_subject_wise       : occupancy rate of each state computed in subject sense
%     model.mean_lifetime_subject_wise                 : mean life time of each state computed in subject sense
%     model.id_of_dominant_states_subject_wise     : id of dominant states subject wise
%     model.id_of_remaining_states                        : id of remaining states for each subject

% taghia@stanford.edu (2016)

if nargin<3
    max_ldim = size(data{1}, 1) - 1;
    opt.n_iter = 100;
    opt.n_init_iter = 10;
    opt.tol = 1e-3;
    opt.noise = 0;
    opt.n_init_learning = 10;
    opt.p_ref = 0.01;
end

if nargin<4
    opt.n_iter = 100;
    opt.n_init_iter = 10;
    opt.tol = 1e-3;
    opt.noise = 0;
    opt.n_init_learning = 10;
    opt.p_ref = 0.01;
end

% modeling
for i = 1: opt.n_init_learning
    display(' ')
    display(['initial learning: ', num2str(i), '/' num2str(opt.n_init_learning)])
    net = vbhafa(data, max_nstates, max_ldim, opt.n_init_iter, opt.tol, opt.noise, 1);
    net_trial{i} = net;
end
for i=1:opt.n_init_learning
    ll(i) = net_trial{i}.Fhist(end,end);
end
display(' ')
display('selecting the best initialization and the best initial model')
[~, id_net] = max(ll);
net_best = net_trial{id_net};
display(' ')
display(['final learning using initial model ', num2str(id_net), '...'])
net_opt = vbhafa(data, max_nstates, max_ldim, opt.n_iter, opt.tol, opt.noise, 1, net_best);
[~, estStatesCell] =  estimateStatesByVitterbi(data, net_opt.params, net_opt.logOutProbs);
[fractional_occupancy_group, mean_life_group]  = compute_occupancy_and_mean_life_group_wise(estStatesCell, max_nstates);
[fractional_occupancy_subj, mean_life_subj]  = compute_occupancy_and_mean_life_subject_wise(estStatesCell, max_nstates);
id_of_dominant_states_group = getDominantStateIdsGroup(estStatesCell, max_nstates);
id_of_dominant_states_subject = getDominantStateIdsSubject(estStatesCell, max_nstates);
display(' ')
display(['final wieghts:', mat2str(sum(net_opt.hidden.Qns))]);
display('done.')

% building model
model.net = net_opt;
model.estimated_covariance = getCovariance(net_opt);
model.estimated_mean = getMean(net_opt);
model.temporal_evolution_of_states = estStatesCell;
model.posterior_probabilities = net_opt.hidden.QnsCell;
model.joint_posterior_probabilities = net_opt.hidden.Qnss;
model.id_of_remaining_states = getRemainingStateIds(model.temporal_evolution_of_states);
model.id_of_dominant_states_group_wise = id_of_dominant_states_group;
model.id_of_dominant_states_subject_wise = id_of_dominant_states_subject;
model.state_transition_probabilities = net_opt.params.stran;
model.fractional_occupancy_group_wise = fractional_occupancy_group;
model.mean_lifetime_group_wise = mean_life_group;
model.fractional_occupancy_subject_wise = fractional_occupancy_subj;
model.mean_lifetime_subject_wise = mean_life_subj;



