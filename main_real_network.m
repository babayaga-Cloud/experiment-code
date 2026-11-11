%% Main Script: Real-world Network Immunization Experiments
clear; clc; close all;
tic;

%% 1. Load Network Data
% Ensure 'example_network.mat' is in the data/ folder
load('Real_Network_SI_data.mat'); 
N = node_info.num_nodes;
E = sum(adj_matrix(:));
out_degrees = sum(adj_matrix, 2);
k = sum(out_degrees) / N;

[~, comp] = graphconncomp(sparse(adj_matrix), 'Directed', true);
scc_sizes = accumarray(comp', 1);
scc_ratio = sum(scc_sizes(scc_sizes > 1)) / N;
max_scc_size = max(scc_sizes) / N;

fprintf('=== Network Basic Information ===\n');
fprintf(' - Nodes: %d\n', N);
fprintf(' - Edges: %d\n', E);
fprintf(' - Average Out-degree: %.1f\n', k);
fprintf(' - SCC Node Ratio: %.2f\n', scc_ratio);
fprintf(' - GSCC Ratio: %.2f\n', max_scc_size);

%% 2. Simulation and Experiment Parameters
beta = 0.25;         % Infection rate
I0 = 0.01;           % Initial infected proportion
Omega = max(1, round(N * I0)); % Initial number of infected nodes

num_trials = 200;                 % Number of independent runs
p_immunize_list = 0:0.1:0.5;      % Immunization proportion array
strategy_names = {'DHI', 'GSCC', 'HDA', 'CI', 'TAD'}; 

% Parameters for Tabu Search (DHI and GSCC strategies)
tabu_params.total_iteration =2000;
tabu_params.num_candidate = 200; 
tabu_params.max_iteration = 200; 
tabu_params.num_tabu = 20; 
ell_CI = 2; % Ball radius for CI strategy

%% 3. Initialize Results Structure
results = struct();
for s_idx = 1:length(strategy_names)
    strat = strategy_names{s_idx};
    results.(strat).mean_infection = [];
    results.(strat).std_infection = [];
    results.(strat).ci95_infection = [];
    results.(strat).scc_ratio = [];
    results.(strat).optimal_strategy = [];
end

%% 4. Main Loop over Immunization Proportions
for idx = 1:length(p_immunize_list)
    p_immunize = p_immunize_list(idx);
    fprintf('\nProcessing immunization proportion: %.0f%%\n', p_immunize * 100);
    
    % 4.1 Find optimal immunization nodes for each strategy
    if p_immunize == 0
        for s_idx = 1:length(strategy_names)
            results.(strategy_names{s_idx}).optimal_strategy = ones(1, N);
        end
    else
        for s_idx = 1:length(strategy_names)
            strategy_name = strategy_names{s_idx};
            switch strategy_name
                case 'HDA'
                    [immunized_nodes, ~] = immunization_strategy('HDA', adj_matrix, p_immunize, []);
                case 'CI'
                    [immunized_nodes, ~] = immunization_strategy('CI', adj_matrix, p_immunize, ell_CI);
                case 'TAD'
                    [immunized_nodes, ~] = immunization_strategy('TAD', adj_matrix, p_immunize);
                case 'GSCC'
                    [immunized_nodes, ~] = immunization_strategy('GSCC', adj_matrix, p_immunize, N, tabu_params);
                case 'DHI'
                    [immunized_nodes, ~] = immunization_strategy('DHI', adj_matrix, p_immunize, Omega, N, tabu_params);
            end
            % Store strategy (0 for immunized, 1 for active)
            optimal_strategy = ones(1, N);
            optimal_strategy(immunized_nodes) = 0;
            results.(strategy_name).optimal_strategy = optimal_strategy;
        end
    end
    
    % 4.2 Run SI epidemic simulation
    fprintf('  Running SI simulations...\n');
    infection_rates = struct();
    for s_idx = 1:length(strategy_names)
        infection_rates.(strategy_names{s_idx}) = zeros(num_trials, 1);
    end
    
    for trial = 1:num_trials
        if p_immunize == 0
            final_infected = directed_SI_simulation(adj_matrix, beta, Omega, N);
            rate = sum(final_infected) / N;
            for s_idx = 1:length(strategy_names)
                infection_rates.(strategy_names{s_idx})(trial) = rate;
            end
        else
            for s_idx = 1:length(strategy_names)
                strat = strategy_names{s_idx};
                immunized_nodes = find(results.(strat).optimal_strategy == 0);
                G_immunized = apply_immunization(adj_matrix, immunized_nodes);
                n_imm = size(G_immunized, 1);
                final_infected = directed_SI_simulation(G_immunized, beta, Omega, n_imm);
                infection_rates.(strat)(trial) = sum(final_infected) / N;
            end
        end
    end
    
    % 4.3 Calculate SCC ratios and store statistics
    for s_idx = 1:length(strategy_names)
        strat = strategy_names{s_idx};
        immunized_nodes = find(results.(strat).optimal_strategy == 0);
        
        % SCC ratio
        scc_ratio = calculate_SCC_Ratio(adj_matrix, immunized_nodes, N);
        results.(strat).scc_ratio = [results.(strat).scc_ratio, scc_ratio];
        
        % Stats
        rates = infection_rates.(strat);
        results.(strat).mean_infection = [results.(strat).mean_infection, mean(rates)];
        results.(strat).std_infection = [results.(strat).std_infection, std(rates)];
        results.(strat).ci95_infection = [results.(strat).ci95_infection, 1.96 * std(rates) / sqrt(num_trials)];
    end
    
    % Print current results
    for s_idx = 1:length(strategy_names)
        strat = strategy_names{s_idx};
        fprintf('  %s: Mean Infection = %.4f, SCC Ratio = %.4f\n', ...
            strat, results.(strat).mean_infection(end), results.(strat).scc_ratio(end));
    end
end

%% 5. Plot Results (User can implement their own custom plotting function here)
% plot_results(results, p_immunize_list, strategy_names, N, k, beta, I0);

%% 6. Save Results
save_filename = 'Real_Network_Immunization_Results.mat';
save(save_filename, 'results', 'p_immunize_list', 'strategy_names');
fprintf('\nResults saved to %s\n', save_filename);
toc;

%% Helper function for SCC Ratio calculation
function scc_ratio = calculate_SCC_Ratio(A, immunized_nodes, N_original)
    if ~isempty(immunized_nodes)
        keep_nodes = setdiff(1:size(A,1), immunized_nodes);
        G_immunized = A(keep_nodes, keep_nodes);
    else
        G_immunized = A;
    end
    if size(G_immunized, 1) == 0
        scc_ratio = 0;
        return;
    end
    [~, comp] = graphconncomp(sparse(G_immunized), 'Directed', true);
    comp_sizes = accumarray(comp', 1);
    significant_sccs = comp_sizes(comp_sizes > 1);
    total_scc_nodes = sum(significant_sccs);
    scc_ratio = total_scc_nodes / N_original;
end
%% Helper: Apply immunization
function G_immunized = apply_immunization(G, immunized_nodes)
    if isempty(immunized_nodes)
        G_immunized = G;
    else
        keep_nodes = setdiff(1:size(G,1), immunized_nodes);
        G_immunized = G(keep_nodes, keep_nodes);
    end
end