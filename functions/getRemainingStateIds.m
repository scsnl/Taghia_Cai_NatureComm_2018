function id_of_remaining_states = getRemainingStateIds(temporal_evolution_of_states)

n_subjs = length(temporal_evolution_of_states);
for subj = 1:n_subjs
      id_of_remaining_states{subj} = unique(temporal_evolution_of_states{subj});
end