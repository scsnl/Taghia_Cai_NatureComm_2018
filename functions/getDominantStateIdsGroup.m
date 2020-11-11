function dominant_states = getDominantStateIdsGroup(temporal_evolution_of_states, max_nstates)
K = max_nstates;
sc =zeros(1,K);
for j=1:K
      sc(j) =  sum(cell2mat(temporal_evolution_of_states)==j);
end
counts_post = sc/sum(sc);
[~, dominant_states] = sort(counts_post,'descend');
