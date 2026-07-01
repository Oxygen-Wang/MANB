function [X_GB, F_GB, history, best_narrow_band] = pso_optimize_narrow_band_full_range(...
    k0, N_layer, P, MaxIter, ...
    k_min_scan, k_max_scan, num_k_coarse, ...
    slope_threshold, width_threshold, amplitude_threshold, ...
    rho, c, d, gamma, xi, epsilon, verbose)
% Optimize narrowest narrow band in 0-10k0 range using PSO algorithm
%
% Parameters:
%   k0: Normalized wavenumber reference (for normalization)
%   N_layer: Sequence length (number of layers)
%   P: Number of particles
%   MaxIter: Maximum number of iterations
%   k_min_scan, k_max_scan: Scanning range (0-10k0)
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
%   F_GB: Global optimal fitness value (narrowest narrow band width)
%   history: Iteration history [iteration number, optimal fitness, maximum difference]
%   best_narrow_band: Optimal narrow band information [k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude]

if nargin < 19
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
best_bands = cell(P, 1);  % Store optimal narrow band information for each particle
parfor p_idx = 1:P
    S_seq = X(p_idx, :);
    [F(p_idx), best_bands{p_idx}] = calc_fitness_full_range(S_seq, k0, k_min_scan, k_max_scan, ...
        num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
end

% Local best
F_LB = F;
X_LB = X;
best_bands_LB = best_bands;

% Global best
[F_GB, idx] = min(F_LB);
X_GB = X_LB(idx, :);
best_narrow_band = best_bands_LB{idx};

if verbose
    fprintf('=== PSO Optimization of Narrowest Narrow Band in 0-10k0 Range ===\n');
    fprintf('k0 = %.6e, Sequence length = %d, Number of particles = %d\n', k0, N_layer, P);
    fprintf('Scanning range: [%.6e, %.6e] (%.2f - %.2f k0)\n', k_min_scan, k_max_scan, k_min_scan/k0, k_max_scan/k0);
    fprintf('Initial optimal narrow band width: %.6e\n', F_GB);
    if ~isempty(best_narrow_band)
        fprintf('Initial optimal narrow band position: k ∈ [%.6e, %.6e] (%.4f - %.4f k0)\n', ...
            best_narrow_band(1), best_narrow_band(2), best_narrow_band(1)/k0, best_narrow_band(2)/k0);
    end
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
    best_bands_new = cell(P, 1);
    parfor k_idx = 1:P
        S_seq = X(k_idx, :);
        [F_new(k_idx), best_bands_new{k_idx}] = calc_fitness_full_range(S_seq, k0, k_min_scan, k_max_scan, ...
            num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d);
    end
    F = F_new;
    
    % Update local best
    for k_idx = 1:P
        if F(k_idx) < F_LB(k_idx)
            F_LB(k_idx) = F(k_idx);
            X_LB(k_idx, :) = X(k_idx, :);
            best_bands_LB{k_idx} = best_bands_new{k_idx};
        end
    end
    
    % Update global best
    [curr_GB, idx] = min(F_LB);
    found_new_best = false;
    if curr_GB < F_GB
        F_GB = curr_GB;
        X_GB = X_LB(idx, :);
        best_narrow_band = best_bands_LB{idx};
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
        if ~isempty(best_narrow_band)
            fprintf('        Narrow band position: k ∈ [%.6e, %.6e] (%.4f - %.4f k0)\n', ...
                best_narrow_band(1), best_narrow_band(2), best_narrow_band(1)/k0, best_narrow_band(2)/k0);
        end
        F_GB_prev = F_GB;
    end
end

history = history(1:iter, :);

if verbose
    fprintf('\n=== Optimization Results ===\n');
    fprintf('Global optimal sequence: ');
    fprintf('%d ', X_GB); fprintf('\n');
    fprintf('Minimum narrow band width Δk/k0 = %.6e\n', F_GB);
    if ~isempty(best_narrow_band)
        fprintf('Optimal narrow band position: k ∈ [%.6e, %.6e] (%.4f - %.4f k0)\n', ...
            best_narrow_band(1), best_narrow_band(2), best_narrow_band(1)/k0, best_narrow_band(2)/k0);
        fprintf('Average slope: %.6e, Maximum slope: %.6e, Amplitude: %.6e\n', ...
            best_narrow_band(4), best_narrow_band(5), best_narrow_band(6));
    end
end

end

%% ================== Subfunction: Calculate Fitness (Find Narrowest Narrow Band in Full Range) ==================
function [fitness, best_band] = calc_fitness_full_range(S_seq, k0, k_min_scan, k_max_scan, ...
    num_k_coarse, slope_threshold, width_threshold, amplitude_threshold, rho, c, d)
% Calculate fitness: find narrowest narrow band in entire scanning range
% If no narrow band is found, return a very large value
% Fitness function: minimize narrow band width

best_band = [];

try
    % Coarse scan entire range
    [real_coarse, ~, k_coarse] = run_scan_from_seq(k_min_scan, k_max_scan, num_k_coarse, S_seq, d, rho, c);
    
    % Detect all narrow bands
    q = real_coarse(1,:);
    dq = diff(q);
    dk = diff(k_coarse);
    slope = abs(dq ./ dk);
    
    % Find high slope regions
    high_slope_mask = slope > slope_threshold;
    high_slope_idx = find(high_slope_mask);
    
    if isempty(high_slope_idx)
        fitness = 1e10;
        return;
    end
    
    % Group contiguous high slope regions
    slope_groups = contiguous_groups(high_slope_idx);
    
    candidates = [];
    for i = 1:size(slope_groups, 1)
        idx_start = max(slope_groups(i, 1), 1);
        idx_end = min(slope_groups(i, 2) + 1, length(k_coarse));
        
        if idx_end > idx_start
            k_min = k_coarse(idx_start);
            k_max = k_coarse(idx_end);
            delta_k_norm = (k_max - k_min) / k0;
            
            % Calculate average slope in this interval
            local_slopes = slope(idx_start:idx_end-1);
            avg_slope = mean(local_slopes);
            max_slope = max(local_slopes);
            
            % Calculate amplitude of band values within narrow band (max - min)
            band_values = q(idx_start:idx_end);
            amplitude = max(band_values) - min(band_values);
            
            % Keep candidates satisfying width and amplitude conditions
            if delta_k_norm < width_threshold && amplitude >= amplitude_threshold
                candidates = [candidates; k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude];
            end
        end
    end
    
    if isempty(candidates)
        fitness = 1e10;
        return;
    end
    
    % Select narrowest narrow band (sort by delta_k_norm)
    [~, sort_idx] = sort(candidates(:,3), 'ascend');
    best_band = candidates(sort_idx(1), :);
    fitness = best_band(3);  % Return normalized narrow band width
    
catch
    % If no narrow band is found, return a very large value
    fitness = 1e10;
    best_band = [];
end

end








