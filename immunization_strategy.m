function [immunized_nodes, additional_output] = immunization_strategy(strategy_name, A, p_immunize, varargin)
% RUN_IMMUNIZATION_STRATEGY Dispatch function to run specified immunization strategy
%
% Inputs:
%   strategy_name - Name of the strategy: 'HDA', 'CI', 'TAD', 'GSCC', 'DHI'
%   A             - Network adjacency matrix
%   p_immunize    - Immunization proportion
%   varargin      - Extra arguments depending on the strategy
%
% Outputs:
%   immunized_nodes   - Indices of immunized nodes
%   additional_output - Extra outputs (varies by strategy)

switch strategy_name
    case 'HDA'
        [immunized_nodes, additional_output] = strategy_HDA(A, p_immunize);
    case 'CI'
        ell = varargin{1};
        [immunized_nodes, additional_output] = strategy_CI(A, ell, p_immunize);
    case 'TAD'
        [immunized_nodes, additional_output] = strategy_TAD(A, p_immunize);
    case 'GSCC'
        n_original = varargin{1};
        tabu_params = varargin{2};
        [strategy_vector, additional_output] = strategy_tabu_GSCC(A, p_immunize, n_original, tabu_params);
        immunized_nodes = find(strategy_vector == 0);
    case 'DHI'
        Omega = varargin{1};
        n_original = varargin{2};
        tabu_params = varargin{3};
        [strategy_vector, additional_output] = strategy_tabu_DHI(A, Omega, p_immunize, n_original, tabu_params);
        immunized_nodes = find(strategy_vector == 0);
    otherwise
        error('Unknown strategy name: %s', strategy_name);
end
end

%% ========================================================================
%  Strategy 1: HDA (High-Degree Adaptive)
%  ========================================================================
function [immunized_nodes, gcc_size_history] = strategy_HDA(A, p_immunize)
    N = size(A, 1);
    num_immunize = round(p_immunize * N);
    immunized_nodes = [];
    gcc_size_history = [];
    original_indices = 1:N; 
    current_A = A;
    
    for step = 1:num_immunize
        current_N = length(original_indices);
        if current_N == 0
            break;
        end
        
        % Calculate sum of in-degree and out-degree
        in_degree = sum(current_A, 1);
        out_degree = sum(current_A, 2)';
        current_degree = in_degree + out_degree;
        
        % Find the node with maximum degree
        [~, idx_in_current] = max(current_degree);
        node_to_remove_original = original_indices(idx_in_current);
        immunized_nodes = [immunized_nodes, node_to_remove_original];
        
        % Update network
        keep_mask = ~ismember(original_indices, immunized_nodes);
        original_indices = original_indices(keep_mask);
        current_A = A(original_indices, original_indices);
        
        % Calculate GSCC size
        if current_N > 1
            [~, comp] = graphconncomp(sparse(current_A), 'Directed', true);
            comp_sizes = accumarray(comp', 1);
            gcc_size = max(comp_sizes);
        else
            gcc_size = 0;
        end
        gcc_size_history = [gcc_size_history, gcc_size];
    end
end

%% ========================================================================
%  Strategy 2: CI (Collective Influence)
%  ========================================================================
function [immunized_nodes, gcc_size_history] = strategy_CI(A, ell, p_immunize)
    N = size(A, 1);
    num_immunize = round(p_immunize * N);
    immunized_nodes = [];
    gcc_size_history = [];
    is_removed = false(N, 1);
    
    for step = 1:num_immunize
        current_nodes = find(~is_removed);
        current_N = length(current_nodes);
        if current_N == 0
            break;
        end
        
        current_A = A(~is_removed, ~is_removed);
        in_degree_current = sum(current_A, 1);
        out_degree_current = sum(current_A, 2)';
        current_degree = in_degree_current + out_degree_current;
        
        ci_values = zeros(current_N, 1);
        for i = 1:current_N
            zi = max(current_degree(i) - 1, 0);
            if zi == 0
                ci_values(i) = 0;
                continue;
            end
            
            u = current_nodes(i);
            frontier_nodes = find_directed_frontier_bfs(A, u, ell, is_removed);
            if isempty(frontier_nodes)
                ci_values(i) = 0;
                continue;
            end
            
            sum_zj = 0;
            for j = 1:length(frontier_nodes)
                v = frontier_nodes(j);
                idx_v = find(current_nodes == v, 1);
                if ~isempty(idx_v)
                    zj = max(current_degree(idx_v) - 1, 0);
                    sum_zj = sum_zj + zj;
                end
            end
            ci_values(i) = zi * sum_zj;
        end
        
        [~, idx_in_current] = max(ci_values);
        node_to_remove = current_nodes(idx_in_current);
        immunized_nodes = [immunized_nodes, node_to_remove];
        is_removed(node_to_remove) = true;
        
        current_A_sub = A(~is_removed, ~is_removed);
        if size(current_A_sub, 1) > 1
            [~, comp] = graphconncomp(sparse(current_A_sub), 'Directed', true);
            comp_sizes = accumarray(comp', 1);
            gcc_size = max(comp_sizes);
        else
            gcc_size = 0;
        end
        gcc_size_history = [gcc_size_history, gcc_size];
    end
end

function frontier_nodes = find_directed_frontier_bfs(A, start_node, ell, is_removed)
    N = size(A, 1);
    visited = false(1, N);
    visited(start_node) = true;
    queue = [start_node];
    dist = zeros(1, N);
    frontier_nodes = [];
    
    while ~isempty(queue)
        v = queue(1);
        queue(1) = [];
        
        if dist(v) >= ell
            if dist(v) == ell
                frontier_nodes = [frontier_nodes, v];
            end
            continue;
        end
        
        neighbors = find(A(v, :) == 1);
        for neighbor = neighbors
            if ~visited(neighbor) && ~is_removed(neighbor)
                visited(neighbor) = true;
                dist(neighbor) = dist(v) + 1;
                queue = [queue, neighbor];
            end
        end
    end
end

%% ========================================================================
%  Strategy 3: TAD (Trophic Analysis-based Dismantling)
%  ========================================================================
function [immunized_nodes, gcc_size_history] = strategy_TAD(A, p_immunize)
    N = size(A, 1);
    num_immunize = round(p_immunize * N);
    immunized_nodes = [];
    gcc_size_history = [];
    is_removed = false(N, 1);
    
    for step = 1:num_immunize
        current_nodes = find(~is_removed);
        current_N = length(current_nodes);
        if current_N == 0
            break;
        end
        
        current_A = A(~is_removed, ~is_removed);
        [~, trophic_levels_current] = compute_trophic_levels(A, is_removed);
        [adjusted_A, ~] = identify_backward_links(current_A, trophic_levels_current);
        NI_centrality = compute_NI_centrality(current_A, adjusted_A, trophic_levels_current);
        
        [~, idx_in_current] = max(NI_centrality);
        node_to_remove = current_nodes(idx_in_current);
        immunized_nodes = [immunized_nodes, node_to_remove];
        is_removed(node_to_remove) = true;
        
        current_A_sub = A(~is_removed, ~is_removed);
        if size(current_A_sub, 1) > 1
            [~, comp] = graphconncomp(sparse(current_A_sub), 'Directed', true);
            comp_sizes = accumarray(comp', 1);
            gcc_size = max(comp_sizes(comp_sizes>1));
        else
            gcc_size = 0;
        end
        gcc_size_history = [gcc_size_history, gcc_size];
    end
end

function [trophic_levels_full, trophic_levels_current] = compute_trophic_levels(A, is_removed)
    N = size(A, 1);
    current_nodes = find(~is_removed);
    current_A = A(current_nodes, current_nodes);
    current_N = length(current_nodes);
    trophic_levels_full = NaN(N, 1);
    
    if current_N == 0
        trophic_levels_current = [];
        return;
    end
    
    in_degree = sum(current_A, 1)';
    out_degree = sum(current_A, 2);
    total_degree = in_degree + out_degree;
    degree_diff = in_degree - out_degree;
    
    u = total_degree;
    A1 = diag(u) - current_A - current_A';
    A1_reg = A1 + 1e-8 * eye(current_N);
    
    h_current = A1_reg \ degree_diff;
    h_current = h_current - min(h_current);
    
    trophic_levels_full(current_nodes) = h_current;
    trophic_levels_current = h_current;
end

function [adjusted_A, backward_links] = identify_backward_links(A, trophic_levels)
    N = size(A, 1);
    adjusted_A = zeros(N, N);
    backward_links = [];
    for i = 1:N
        for j = 1:N
            if A(i, j) == 1
                if trophic_levels(j) - trophic_levels(i) < 0
                    adjusted_A(i, j) = 1;
                    adjusted_A(j, i) = 1;
                    backward_links = [backward_links; i, j];
                end
            end
        end
    end
end

function NI_centrality = compute_NI_centrality(~, adjusted_A, trophic_levels)
    N = size(adjusted_A, 1);
    NI_centrality = zeros(N, 1);
    for i = 1:N
        sum_value = 0;
        neighbors = find(adjusted_A(i, :) == 1);
        for j = neighbors
            weight = (trophic_levels(j) - trophic_levels(i) - 1)^2;
            sum_value = sum_value + weight;
        end
        NI_centrality(i) = sum_value;
    end
end

%% ========================================================================
%  Strategy 4: Directed Tabu Search (DTS) for GSCC Minimization
%  ========================================================================
function [optimal_strategy, optimal_value] = strategy_tabu_GSCC(G, p_immunize, n_original, tabu_params)
    n = size(G, 1);
    num_immunize = round(p_immunize * n);
    total_iteration = tabu_params.total_iteration;
    num_candidate = tabu_params.num_candidate;
    max_iteration = tabu_params.max_iteration;
    num_tabu = tabu_params.num_tabu;
    
    % ---------- Initial solution ----------
    [scc_ids_orig, scc_ri_orig, sorted_scc_orig] = compute_scc_influence(G, n_original);
    num_scc_orig = length(scc_ri_orig);
    scc_members_orig = cell(num_scc_orig, 1);
    for i = 1:n
        scc_members_orig{scc_ids_orig(i)} = [scc_members_orig{scc_ids_orig(i)}, i];
    end
    
    x_current = ones(1, n);
    selected_count = 0;
    for idx = 1:num_scc_orig
        scc = sorted_scc_orig(idx);
        members = scc_members_orig{scc};
        perm = randperm(length(members));
        for m = 1:length(members)
            if selected_count >= num_immunize, break; end
            node = members(perm(m));
            if x_current(node) == 1
                x_current(node) = 0;
                selected_count = selected_count + 1;
            end
        end
        if selected_count >= num_immunize, break; end
    end
    while selected_count < num_immunize
        unselected = find(x_current == 1);
        if isempty(unselected), break; end
        x_current(unselected(randi(length(unselected)))) = 0;
        selected_count = selected_count + 1;
    end
    
    y_current = calculate_GSCC_for_strategy(G, n_original, x_current);
    x_optimal = x_current; y_optimal = y_current;
    Tlist = []; num_iteration = 0; total_i = 0;
    best_x = x_optimal; best_y = y_optimal;
    
    % ---------- Main Loop ----------
    while (num_iteration < max_iteration) && (total_i < total_iteration)
        total_i = total_i + 1;
        
        current_nodes = find(x_current == 1);
        if isempty(current_nodes), break; end
        current_A = G(current_nodes, current_nodes);
        [scc_ids_local, scc_ri_local, sorted_scc_local] = compute_scc_influence(current_A, n_original);
        num_scc_local = length(scc_ri_local);
        
        scc_members = cell(num_scc_local, 1);
        for i = 1:length(current_nodes)
            node_orig = current_nodes(i);
            scc_members{scc_ids_local(i)} = [scc_members{scc_ids_local(i)}, node_orig];
        end
        
        node_to_sub = zeros(n, 1);
        node_to_sub(current_nodes) = 1:length(current_nodes);
        
        tabu_temp = zeros(1, num_candidate);
        x_candidate = zeros(num_candidate, n);
        y_candidate = zeros(1, num_candidate);
        
        for k = 1:num_candidate
            x_temp = x_current;
            
            scc_idx = randi(num_scc_local);
            scc = sorted_scc_local(scc_idx);
            members = scc_members{scc};
            if isempty(members), continue; end
            start_node = members(randi(length(members)));
            
            steps = randi([1, 3]);
            v = start_node;
            for step = 1:steps
                sub_v = node_to_sub(v);
                if sub_v == 0, break; end
                out_neighbors_sub = find(current_A(sub_v, :) == 1);
                if isempty(out_neighbors_sub), break; end
                out_neighbors = current_nodes(out_neighbors_sub);
                v = out_neighbors(randi(length(out_neighbors)));
            end
            
            if x_temp(v) == 0
                non_immunized = find(x_temp == 1);
                if isempty(non_immunized), continue; end
                u = non_immunized(randi(length(non_immunized)));
                x_temp(v) = 1; x_temp(u) = 0;
            else
                immunized = find(x_temp == 0);
                if isempty(immunized), continue; end
                u = immunized(randi(length(immunized)));
                x_temp(v) = 0; x_temp(u) = 1;
            end
            
            x_candidate(k, :) = x_temp;
            tabu_temp(k) = iftabu(Tlist, x_temp);
            y_candidate(k) = calculate_GSCC_for_strategy(G, n_original, x_temp);
        end
        
        [y_min, num_min] = min(y_candidate);
        if y_min < y_optimal
            x_current = x_candidate(num_min, :); y_current = y_min;
            x_optimal = x_current; y_optimal = y_current;
            num_iteration = 0;
            if y_optimal < best_y
                best_x = x_optimal; best_y = y_optimal;
            end
        else
            candidate_indices = find(tabu_temp == 1);
            if ~isempty(candidate_indices)
                [y_min2, min2_idx] = min(y_candidate(candidate_indices));
                x_current = x_candidate(candidate_indices(min2_idx), :); y_current = y_min2;
            else
                [y_min2, min2_idx] = min(y_candidate);
                x_current = x_candidate(min2_idx, :); y_current = y_min2;
            end
            num_iteration = num_iteration + 1;
        end
        
        if size(Tlist, 1) < num_tabu
            Tlist = [Tlist; x_current];
        else
            Tlist(1, :) = []; Tlist = [Tlist; x_current];
        end
    end
    optimal_strategy = best_x; optimal_value = best_y;
end

%% ========================================================================
%  Strategy 5: Directed Tabu Search (DTS) for DHI Minimization
%  ========================================================================
function [optimal_strategy, optimal_value] = strategy_tabu_DHI(G, Omega, p_immunize, n_original, tabu_params)
    n = size(G, 1);
    num_immunize = round(p_immunize * n);
    total_iteration = tabu_params.total_iteration;
    num_candidate = tabu_params.num_candidate;
    max_iteration = tabu_params.max_iteration;
    num_tabu = tabu_params.num_tabu;
    
    % ---------- Initial solution ----------
    [scc_ids_orig, scc_ri_orig, sorted_scc_orig] = compute_scc_influence(G, n_original);
    num_scc_orig = length(scc_ri_orig);
    scc_members_orig = cell(num_scc_orig, 1);
    for i = 1:n
        scc_members_orig{scc_ids_orig(i)} = [scc_members_orig{scc_ids_orig(i)}, i];
    end
    
    x_current = ones(1, n);
    selected_count = 0;
    for idx = 1:num_scc_orig
        scc = sorted_scc_orig(idx);
        members = scc_members_orig{scc};
        perm = randperm(length(members));
        for m = 1:length(members)
            if selected_count >= num_immunize, break; end
            node = members(perm(m));
            if x_current(node) == 1
                x_current(node) = 0;
                selected_count = selected_count + 1;
            end
        end
        if selected_count >= num_immunize, break; end
    end
    while selected_count < num_immunize
        unselected = find(x_current == 1);
        if isempty(unselected), break; end
        x_current(unselected(randi(length(unselected)))) = 0;
        selected_count = selected_count + 1;
    end
    
    y_current = calculate_DHI_for_strategy(G, Omega, n_original, x_current);
    x_optimal = x_current; y_optimal = y_current;
    Tlist = []; num_iteration = 0; total_i = 0;
    best_x = x_optimal; best_y = y_optimal;
    
    % ---------- Main Loop ----------
    while (num_iteration < max_iteration) && (total_i < total_iteration)
        total_i = total_i + 1;
        
        current_nodes = find(x_current == 1);
        if isempty(current_nodes), break; end
        current_A = G(current_nodes, current_nodes);
        [scc_ids_local, scc_ri_local, sorted_scc_local] = compute_scc_influence(current_A, n_original);
        num_scc_local = length(scc_ri_local);
        
        scc_members = cell(num_scc_local, 1);
        for i = 1:length(current_nodes)
            node_orig = current_nodes(i);
            scc_members{scc_ids_local(i)} = [scc_members{scc_ids_local(i)}, node_orig];
        end
        
        node_to_sub = zeros(n, 1);
        node_to_sub(current_nodes) = 1:length(current_nodes);
        
        tabu_temp = zeros(1, num_candidate);
        x_candidate = zeros(num_candidate, n);
        y_candidate = zeros(1, num_candidate);
        
        for k = 1:num_candidate
            x_temp = x_current;
            
            scc_idx = randi(num_scc_local);
            scc = sorted_scc_local(scc_idx);
            members = scc_members{scc};
            if isempty(members), continue; end
            start_node = members(randi(length(members)));
            
            steps = randi([1, 3]);
            v = start_node;
            for step = 1:steps
                sub_v = node_to_sub(v);
                if sub_v == 0, break; end
                out_neighbors_sub = find(current_A(sub_v, :) == 1);
                if isempty(out_neighbors_sub), break; end
                out_neighbors = current_nodes(out_neighbors_sub);
                v = out_neighbors(randi(length(out_neighbors)));
            end
            
            if x_temp(v) == 0
                non_immunized = find(x_temp == 1);
                if isempty(non_immunized), continue; end
                u = non_immunized(randi(length(non_immunized)));
                x_temp(v) = 1; x_temp(u) = 0;
            else
                immunized = find(x_temp == 0);
                if isempty(immunized), continue; end
                u = immunized(randi(length(immunized)));
                x_temp(v) = 0; x_temp(u) = 1;
            end
            
            x_candidate(k, :) = x_temp;
            tabu_temp(k) = iftabu(Tlist, x_temp);
            y_candidate(k) = calculate_DHI_for_strategy(G, Omega, n_original, x_temp);
        end
        
        [y_min, num_min] = min(y_candidate);
        if y_min < y_optimal
            x_current = x_candidate(num_min, :); y_current = y_min;
            x_optimal = x_current; y_optimal = y_current;
            num_iteration = 0;
            if y_optimal < best_y
                best_x = x_optimal; best_y = y_optimal;
            end
        else
            candidate_indices = find(tabu_temp == 1);
            if ~isempty(candidate_indices)
                [y_min2, min2_idx] = min(y_candidate(candidate_indices));
                x_current = x_candidate(candidate_indices(min2_idx), :); y_current = y_min2;
            else
                [y_min2, min2_idx] = min(y_candidate);
                x_current = x_candidate(min2_idx, :); y_current = y_min2;
            end
            num_iteration = num_iteration + 1;
        end
        
        if size(Tlist, 1) < num_tabu
            Tlist = [Tlist; x_current];
        else
            Tlist(1, :) = []; Tlist = [Tlist; x_current];
        end
    end
    optimal_strategy = best_x; optimal_value = best_y;
end

%% ========================================================================
%  Helper Functions (SCC and Metrics)
%  ========================================================================

function [scc_ids, scc_ri, sorted_scc_indices] = compute_scc_influence(G, n_original)
    n = size(G,1);
    if n == 0
        scc_ids = []; scc_ri = []; sorted_scc_indices = [];
        return;
    end
    G_di = digraph(G);
    scc_ids = conncomp(G_di, 'Type', 'strong');
    num_scc = max(scc_ids);
    scc_sizes = accumarray(scc_ids', 1);
    
    [u, v] = find(G);
    src_scc = scc_ids(u);
    dst_scc = scc_ids(v);
    mask = src_scc ~= dst_scc;
    cond_adj = sparse(src_scc(mask), dst_scc(mask), 1, num_scc, num_scc);
    cond_G = digraph(cond_adj);
    
    reachable_counts = zeros(num_scc, 1);
    if num_scc > 1
        dist_matrix = distances(cond_G);
        for i = 1:num_scc
            reachable_scc = find(dist_matrix(i,:) < Inf);
            reachable_counts(i) = sum(scc_sizes(reachable_scc));
        end
    else
        reachable_counts = scc_sizes;
    end
    scc_ri = reachable_counts / n_original;
    [~, sorted_indices] = sort(scc_ri, 'descend');
    sorted_scc_indices = sorted_indices;
end

function value = calculate_GSCC_for_strategy(G, n_original, strategy)
    immunized_nodes = find(strategy == 0);
    if isempty(immunized_nodes)
        G_temp = G;
    else
        keep = setdiff(1:size(G,1), immunized_nodes);
        G_temp = G(keep, keep);
    end
    value = calculate_GSCC(G_temp, n_original);
end

function value = calculate_DHI_for_strategy(G, Omega, n_original, strategy)
    immunized_nodes = find(strategy == 0);
    if isempty(immunized_nodes)
        G_temp = G;
    else
        keep = setdiff(1:size(G,1), immunized_nodes);
        G_temp = G(keep, keep);
    end
    value = calculate_DHI(G_temp, Omega, n_original);
end

function gcc_ratio = calculate_GSCC(G, n_original)
    if ~isa(G, 'digraph')
        G_di = digraph(G);
    else
        G_di = G;
    end
    n = size(G, 1);
    if n == 0
        gcc_ratio = 0; return; 
    end
    bin = conncomp(G_di, 'Type', 'strong');
    comp_sizes = accumarray(bin', 1);
    if isempty(comp_sizes(comp_sizes>1))
        gcc_ratio = 0;
    else
        gcc_ratio = max(comp_sizes) / n_original;
    end
end

function value = calculate_DHI(G_temp, Omega, n_original)
    n_curr = size(G_temp,1);
    if n_curr == 0
        value = 0; return;
    end
    G_di = digraph(G_temp);
    bin = conncomp(G_di, 'Type', 'strong');
    num_sccs = max(bin);
    scc_sizes = accumarray(bin', 1);
    
    edge_list = table2array(G_di.Edges);
    u_ids = edge_list(:, 1);
    v_ids = edge_list(:, 2);
    src_scc = bin(u_ids);
    dst_scc = bin(v_ids);
    
    mask = src_scc ~= dst_scc;
    cond_adj = sparse(src_scc(mask), dst_scc(mask), 1, num_sccs, num_sccs);
    cond_G = digraph(cond_adj);
    dist_matrix = distances(cond_G); 
    
    reachable_sizes = zeros(num_sccs, 1);
    for i = 1:num_sccs
        reachable_scc_ids = find( dist_matrix(i, :) < Inf);
        reachable_sizes(i) = sum(scc_sizes(reachable_scc_ids));
    end
    
    if isempty(scc_sizes(scc_sizes>1))
        value = 0;
    else
        p_i = scc_sizes / n_original;          
        r_i = reachable_sizes / n_original;    
        value = sum(p_i .* (1 - (1 - p_i).^(Omega)) .* r_i);
    end
end

function result = iftabu(Tlist, x_candidate)
    if isempty(Tlist)
        result = 1; return;
    end
    for i = 1:size(Tlist, 1)
        if isequal(Tlist(i, :), x_candidate)
            result = 0; return;
        end
    end
    result = 1;
end