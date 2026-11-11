function final_infected = fast_directed_SI_simulation(G, beta, Omega, n)
    % FAST_DIRECTED_SI_SIMULATION Fast simulation of the SI model on a directed network
    % Inputs:
    %   G     - Adjacency matrix
    %   beta  - Infection rate
    %   Omega - Initial number of infected nodes
    %   n     - Current number of nodes in the network
    % Output:
    %   final_infected - Boolean vector indicating final infection status

    num_initial = min(Omega, n);
    if num_initial == 0
        final_infected = false(n, 1);
        return;
    end
    
    % Initialize infection states
    infected = zeros(n, 1);
    initial_infected = randperm(n, num_initial);
    infected(initial_infected) = true;
    new_infected = infected;
    susceptible = ~infected;
    
    % Pre-compute sparse adjacency matrix for fast matrix multiplication
    G_sparse = sparse(G > 0);
    
    for step = 1:1000
        % Calculate infection pressure from all in-neighbors
        infection_pressure = (G_sparse)' * double(new_infected);
        infection_prob = 1 - (1 - beta).^infection_pressure;
        
        % Determine new infections
        new_infections = (rand(n, 1) < infection_prob) & susceptible;
        
        if ~any(new_infections)
            break;
        end
        
        % Update states
        new_infected = new_infected | new_infections;
        susceptible = ~new_infected;
    end
    final_infected = new_infected;
end