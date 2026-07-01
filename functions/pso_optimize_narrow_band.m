function [X_GB, F_GB, history] = pso_optimize_narrow_band(k0, N_layer, P, MaxIter, ...
    k_min0, k_max0, num_k_coarse, slope_threshold, width_threshold, ...
    amplitude_threshold, rho, c, d, gamma, xi, epsilon, verbose)
% Optimize narrow band width at k0 to minimum using PSO algorithm
%
% Parameters:
%   k0: Normalized wavenumber reference
%   N_layer: Sequence length (number of layers)
%   P: Number of particles
%   MaxIter: Maximum number of iterations
%   k_min0, k_max0: Scanning range
%   num_k_coarse: Number of coarse scanning points
%   slope_threshold: Slope threshold
%   width_threshold: Width threshold
%   amplitude_threshold: Amplitude threshold
%   rho: Density
%   c: Sound speed
%   d: Total length
%   gamma: Inertia weight adjustment parameter
%   xi: Small positive number
%   epsilon: Convergence threshold
%   verbose: Whether to display detailed information (default true)
%
% Returns:
%   X_GB: Global optimal sequence
%   F_GB: Global optimal fitness value (narrow band width)
%   history: Iteration history [iteration number, optimal fitness, maximum difference]

if nargin < 18
    verbose = true;
end

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir));

%% ====== Initialize Particles ======
% Particle positions: cross-sectional area per layer can only be 1 or 2
X = randi([1, 2], P, N_layer);

% Particle velocities (can be continuous)
V = zeros(P, N_layer);

% Calculate initial fitness (narrow band width) - parallelized using parfor
F = zeros(P, 1);
parfor p_idx = 1:P
    S_seq = X(p_idx, :);
    F(p_idx) = calc_fitness(S_seq, k0, k_min0, k_max0, num_k_coarse, ...
        slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
end

% Local best
F_LB = F;
X_LB = X;

% Global best
[F_GB, idx] = min(F_LB);
X_GB = X_LB(idx, :);

if verbose
    fprintf('=== PSO Optimization of Narrow Band Width at k0 ===\n');
    fprintf('k0 = %.6e, Sequence length = %d, Number of particles = %d\n', k0, N_layer, P);
    fprintf('Initial optimal narrow band width: %.6e\n', F_GB);
end

% Record previous optimal value to determine if new optimal solution is found
F_GB_prev = F_GB;

%% ====== Iterative Optimization ======
converged = false;
iter = 0;
history = zeros(MaxIter, 3); % [iteration number, optimal fitness, maximum difference]

while iter < MaxIter && ~converged
    iter = iter + 1;
    
    % Pre-calculate velocity and position updates for all particles (serial part)
    for k_idx = 1:P
        % ---------- Adaptive parameter calculation ----------
        c1_k = exp(F_GB) / (exp(F_GB) + exp(F_LB(k_idx)));
        c2_k = exp(F_LB(k_idx)) / (exp(F_GB) + exp(F_LB(k_idx)));
        
        if F_GB < F(k_idx) && F(k_idx) < F_LB(k_idx)
            omega_k = xi * (F(k_idx) / F_GB)^gamma;
        else
            omega_k = (2*F(k_idx) - F_GB - F_LB(k_idx))^gamma;
        end
        
        % Ensure omega_k is within reasonable range
        omega_k = max(0.1, min(0.9, omega_k));
        
        % ---------- Velocity update ----------
        r1 = rand(1, N_layer);
        r2 = rand(1, N_layer);
        V(k_idx, :) = omega_k * V(k_idx, :) + c1_k * r1 .* (X_LB(k_idx, :) - X(k_idx, :)) ...
                                   + c2_k * r2 .* (X_GB - X(k_idx, :));
        
        % Limit velocity range
        V_max = 2.0;
        V(k_idx, :) = max(-V_max, min(V_max, V(k_idx, :)));
        
        % ---------- Position update ----------
        X(k_idx, :) = X(k_idx, :) + V(k_idx, :);
        
        % ---------- Discretize to 1 or 2 ----------
        X(k_idx, :) = round(X(k_idx, :));
        X(k_idx, :) = max(min(X(k_idx, :), 2), 1);
    end
    
    % Parallel calculate fitness for all particles
    F_new = zeros(P, 1);
    parfor k_idx = 1:P
        S_seq = X(k_idx, :);
        F_new(k_idx) = calc_fitness(S_seq, k0, k_min0, k_max0, num_k_coarse, ...
            slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
    end
    F = F_new;
    
    % Update local best
    for k_idx = 1:P
        if F(k_idx) < F_LB(k_idx)
            F_LB(k_idx) = F(k_idx);
            X_LB(k_idx, :) = X(k_idx, :);
        end
    end
    
    % Update global best
    [curr_GB, idx] = min(F_LB);
    found_new_best = false;
    if curr_GB < F_GB
        F_GB = curr_GB;
        X_GB = X_LB(idx, :);
        found_new_best = true;
    end
    
    % Convergence check
    max_diff = max(abs(F_LB - F_GB));
    history(iter, :) = [iter, F_GB, max_diff];
    
    if max_diff < epsilon
        converged = true;
        if verbose
            fprintf('Algorithm converged after %d iterations\n', iter);
        end
    end
    
    % Only output when new optimal solution is found
    if verbose && found_new_best
        fprintf('Iteration %3d: Found new optimal solution, narrow band width = %.6e (improvement: %.6e)\n', ...
                iter, F_GB, F_GB_prev - F_GB);
        F_GB_prev = F_GB;
    end
end

history = history(1:iter, :);

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

