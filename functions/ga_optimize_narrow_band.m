function [X_GB, F_GB, history] = ga_optimize_narrow_band(k0, N_layer, pop_size, MaxGen, ...
    k_min0, k_max0, num_k_coarse, slope_threshold, width_threshold, ...
    amplitude_threshold, rho, c, d, pc, pm, elite_ratio, verbose)
% Optimize narrow band width at k0 to minimum using Genetic Algorithm (GA)
%
% Parameters:
%   k0: Normalized wavenumber reference
%   N_layer: Sequence length (number of layers)
%   pop_size: Population size
%   MaxGen: Maximum number of generations
%   k_min0, k_max0: Scanning range
%   num_k_coarse: Number of coarse scanning points
%   slope_threshold: Slope threshold
%   width_threshold: Width threshold
%   amplitude_threshold: Amplitude threshold
%   rho: Density
%   c: Sound speed
%   d: Total length
%   pc: Crossover probability (default 0.8)
%   pm: Mutation probability (default 0.1)
%   elite_ratio: Elite retention ratio (default 0.1)
%   verbose: Whether to display detailed information (default true)
%
% Returns:
%   X_GB: Global optimal sequence
%   F_GB: Global optimal fitness value (narrow band width)
%   history: Iteration history [generation, optimal fitness, average fitness, maximum difference]

if nargin < 13
    pc = 0.8;  % Crossover probability
end
if nargin < 14
    pm = 0.1;  % Mutation probability
end
if nargin < 15
    elite_ratio = 0.1;  % Elite retention ratio
end
if nargin < 16
    verbose = true;
end

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir));

%% ====== Initialize Population ======
% Population: each individual is a sequence, cross-sectional area per layer can only be 1 or 2
pop = randi([1, 2], pop_size, N_layer);

% Calculate initial fitness (narrow band width)
fitness = zeros(pop_size, 1);
for i = 1:pop_size
    S_seq = pop(i, :);
    fitness(i) = calc_fitness(S_seq, k0, k_min0, k_max0, num_k_coarse, ...
        slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
end

% Global best
[F_GB, best_idx] = min(fitness);
X_GB = pop(best_idx, :);

if verbose
    fprintf('=== GA Optimization of Narrow Band Width at k0 ===\n');
    fprintf('k0 = %.6e, Sequence length = %d, Population size = %d\n', k0, N_layer, pop_size);
    fprintf('Initial optimal narrow band width: %.6e\n', F_GB);
end

%% ====== Iterative Optimization ======
history = zeros(MaxGen, 4); % [generation, optimal fitness, average fitness, maximum difference]
elite_size = max(1, round(pop_size * elite_ratio));  % Number of elite individuals

for gen = 1:MaxGen
    % Record current generation information
    avg_fitness = mean(fitness);
    max_fitness = max(fitness);
    min_fitness = min(fitness);
    max_diff = max_fitness - min_fitness;
    history(gen, :) = [gen, min_fitness, avg_fitness, max_diff];
    
    % Update global best
    [curr_min, curr_best_idx] = min(fitness);
    if curr_min < F_GB
        F_GB = curr_min;
        X_GB = pop(curr_best_idx, :);
    end
    
    if verbose && (mod(gen, 10) == 0 || gen == 1 || gen == MaxGen)
        fprintf('Generation %3d: Optimal narrow band width = %.6e, Average = %.6e, Maximum difference = %.8e\n', ...
                gen, F_GB, avg_fitness, max_diff);
    end
    
    % Convergence check
    if max_diff < 1e-6
        if verbose
            fprintf('Algorithm converged after generation %d\n', gen);
        end
        break;
    end
    
    % ====== Selection Operation (Roulette Wheel Selection) ======
    % Convert fitness to selection probability (smaller fitness is better, so take reciprocal)
    fitness_inv = 1 ./ (fitness + 1e-10);  % Avoid division by zero
    prob = fitness_inv / sum(fitness_inv);
    cum_prob = cumsum(prob);
    
    % Elite retention
    [~, elite_idx] = sort(fitness, 'ascend');
    elite_pop = pop(elite_idx(1:elite_size), :);
    elite_fitness = fitness(elite_idx(1:elite_size));
    
    % Select new population (retain elite)
    new_pop = zeros(pop_size, N_layer);
    new_fitness = zeros(pop_size, 1);
    new_pop(1:elite_size, :) = elite_pop;
    new_fitness(1:elite_size) = elite_fitness;
    
    % Roulette wheel selection for remaining individuals
    for i = (elite_size + 1):pop_size
        r = rand;
        selected_idx = find(cum_prob >= r, 1);
        if isempty(selected_idx)
            selected_idx = pop_size;
        end
        new_pop(i, :) = pop(selected_idx, :);
        new_fitness(i) = fitness(selected_idx);
    end
    
    % ====== Crossover Operation ======
    % Random pairing (elite do not participate in crossover)
    pair_indices = randperm(pop_size - elite_size) + elite_size;
    pair_indices = reshape(pair_indices(1:2*floor((pop_size-elite_size)/2)), [], 2);
    
    for p = 1:size(pair_indices, 1)
        if rand < pc  % Crossover probability
            idx1 = pair_indices(p, 1);
            idx2 = pair_indices(p, 2);
            
            % Single-point crossover
            cross_point = randi([1, N_layer-1]);
            temp = new_pop(idx1, cross_point+1:end);
            new_pop(idx1, cross_point+1:end) = new_pop(idx2, cross_point+1:end);
            new_pop(idx2, cross_point+1:end) = temp;
            
            % Recalculate fitness
            new_fitness(idx1) = calc_fitness(new_pop(idx1, :), k0, k_min0, k_max0, ...
                num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
            new_fitness(idx2) = calc_fitness(new_pop(idx2, :), k0, k_min0, k_max0, ...
                num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
        end
    end
    
    % ====== Mutation Operation ======
    for i = (elite_size + 1):pop_size
        if rand < pm  % Mutation probability
            % Randomly select mutation position
            mut_pos = randi(N_layer);
            % Flip value at this position (1 becomes 2, 2 becomes 1)
            new_pop(i, mut_pos) = 3 - new_pop(i, mut_pos);
            % Recalculate fitness
            new_fitness(i) = calc_fitness(new_pop(i, :), k0, k_min0, k_max0, ...
                num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
        end
    end
    
    % Update population
    pop = new_pop;
    fitness = new_fitness;
end

history = history(1:gen, :);

if verbose
    fprintf('\n=== Optimization Results ===\n');
    fprintf('Global optimal sequence: ');
    fprintf('%d ', X_GB); fprintf('\n');
    fprintf('Minimum narrow band width Δk/k0 = %.6e\n', F_GB);
end

end

%% ================== Subfunction: Calculate Fitness (Narrow Band Width) ==================
function fitness = calc_fitness(S_seq, k0, k_min0, k_max0, num_k_coarse, ...
    slope_threshold, width_threshold, amplitude_threshold, rho, c, d)
% Calculate fitness: find narrow band width near k0
% If no narrow band is found, return a very large value
% Fitness function: minimize narrow band width, while considering whether narrow band is near k0

try
    % Use find_narrow_band_from_seq to find narrow band
    [final_result, ~, candidates, ~, ~] = find_narrow_band_from_seq(...
        S_seq, d, k0, ...
        k_min0, k_max0, ...
        num_k_coarse, 100, ...  % Number of coarse scanning points, number of fine scanning points
        slope_threshold, width_threshold, amplitude_threshold, ...
        20, 0.3, 1.5, ...  % max_refinement_iter, refinement_margin, refinement_points_factor
        rho, c, ...
        false);  % verbose = false
    
    % If there are multiple candidate narrow bands, select the one closest to k0
    if ~isempty(candidates) && size(candidates, 1) > 1
        % Calculate center frequency of each candidate narrow band
        k_centers = (candidates(:, 1) + candidates(:, 2)) / 2;
        % Find narrow band closest to k0
        [~, closest_idx] = min(abs(k_centers - k0));
        delta_k_norm = candidates(closest_idx, 3);
        k_center = k_centers(closest_idx);
    else
        % Use final result
        delta_k_norm = final_result(3);
        k_center = (final_result(1) + final_result(2)) / 2;
    end
    
    % Fitness is normalized narrow band width Δk/k0
    % If narrow band center is far from k0, add penalty term
    distance_penalty = abs(k_center - k0) / k0;
    if distance_penalty > 0.1  % If distance from k0 exceeds 10%, add penalty
        fitness = delta_k_norm * (1 + 10 * distance_penalty);
    else
        fitness = delta_k_norm;
    end
    
catch
    % If no narrow band is found, return a very large value
    fitness = 1e10;
end

end








