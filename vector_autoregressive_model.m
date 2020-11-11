function [stats_group, stats_subj] = vector_autoregressive_model(data, group_model, subject_model, n_iter, p_ref_value)
%Inputs: 
%	data 				       : {data_subj1, data_subj2,...,data_subjS}, data_subj1 is a D-by-N matrix where D is dimension 
%     group_model         : BSDS group_model
%    subject_model        :  model for each subject which is computed from  group model (optional)
%     n_iter                   : number of iterations
%     p_value                : p_value

%outputs:
%    stats.Pvals             : p_vlaues
%    stats.Acoeffs          : AR coeffs
%    stats.Zscores          : z-scores
%    stats.Map               : MAP

if nargin<3
      n_iter=100;
      p_ref_value = 0.01;
      subject_model = [];
elseif nargin<4
      p_ref_value = 0.01;
      n_iter=100;
elseif nargin<5
      p_ref_value = 0.01;
else
end
      
stats_group =  posthocVARfromData(data, group_model.net, cell2mat(group_model.temporal_evolution_of_states), 0, n_iter, p_ref_value);
if isempty(subject_model)==0
      n_subjs = length(group_model.temporal_evolution_of_states);
      for subj=1:n_subjs
            stats_subj{subj} =  posthocVARfromData(data(subj), subject_model.net{subj},  group_model.temporal_evolution_of_states{subj}, 0, n_iter, p_ref_value);
      end
else
      stats_subj = [];
end
