
% this script requires posterior probabilities learned from BSDSS model
%Inputs:
% group_model                            : trained model in group level
% data                                        : group data as a cell of individual subject matrices D-by-N
% maxIter(optional)                     : number of iterations (default 60)
% pca_flag                                   : if group_net contains information about pca_flag this input should be discarded however, if the the pca_flag is not saved in the trained model it should be provided as either 0 or 1. pca_flag=0 means the noise variance is allowed to vary across dimensions and pca_flag=1 means the noise variance is the same across dimensions
%Outputs:
% model.net                                 : subject level BSDSS net
% model.stats                               : subject level AR stats (Z-scores and p-values)
% model.estimated_covariance       : subject level estimated covariance matrices for each state

% taghia@stanford.edu (2016)


function model = compute_subject_level_stats(data, group_model, maxIter, pca_flag)
if nargin<3
      maxIter = 60;
      pca_flag = group_model.net.hparams.pcaflag;
end
if nargin<4
      pca_flag = group_model.net.hparams.pcaflag;
end

nSubjs = length(data);
maxLocalDim = size(group_model.net.params.Lm{1}, 2)-1;
for subj = 1:nSubjs
      display(' ');
      display(['subject :', num2str(subj)]);
      net_subj = vbhafa_z(data(subj), group_model.net.hidden.QnsCell{subj}, maxLocalDim, maxIter, 1e-3 , pca_flag);
      estCov{subj} = getCovariance(net_subj);
      estMean{subj} = getMean(net_subj);
      net{subj} = net_subj;
end
display(' ');
model.net = net;
model.estimated_covarinace = estCov;
model.estimated_mean = estMean;
