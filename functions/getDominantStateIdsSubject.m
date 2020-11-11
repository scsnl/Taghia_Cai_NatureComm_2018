function dominant_states = getDominantStateIdsSubject(temporal_evolution_of_states, max_nstates)
K = max_nstates;
n_subjs = length(temporal_evolution_of_states);
for subj=1:n_subjs
      sc =zeros(1,K);
      for j=1:K
            sc(j) =  sum(temporal_evolution_of_states{subj}==j);
      end
      counts_post = sc/sum(sc);
      [~, dominant_states(subj, :)] = sort(counts_post,'descend');
end
