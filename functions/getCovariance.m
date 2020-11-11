function cov_state = getCovariance(net)
nStates = length(net.params.Lm);
for state=1:nStates
      cov_state{state}= net.params.Lm{state}(:,2:end)*net.params.Lm{state}(:,2:end)' + diag(net.hparams.psii)^-1;
end
