function mean_state = getMean(net)
nStates = length(net.params.Lm);
for state=1:nStates
      mean_state{state}= net.params.Lm{state}(:,1);
end
