function [fractional_occupancy, mean_life]  = compute_occupancy_and_mean_life_group_wise(temporal_evolution_of_states, max_nstates)

[fractional_occupancy, mean_life,~] = summary_stats_fast(cell2mat(temporal_evolution_of_states), 1:max_nstates);

