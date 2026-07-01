function [final_result, refinement_history, candidates, k_coarse, real_coarse] = find_narrow_band_from_seq(...
    S_seq, d, k0, ...
    k_min0, k_max0, ...
    num_k_coarse, num_k_fine, ...
    slope_threshold, width_threshold, amplitude_threshold, ...
    max_refinement_iter, refinement_margin, refinement_points_factor, ...
    rho, c, ...
    verbose)
% Main function to find narrow bands from sequence array
%
% Parameters:
%   S_seq: Cross-section sequence array
%   d: Total length
%   k0: Normalized wavenumber reference
%   k_min0, k_max0: Scanning range
%   num_k_coarse: Number of coarse scanning points
%   num_k_fine: Number of fine scanning points
%   slope_threshold: Absolute value threshold for slope
%   width_threshold: Δk/k0 threshold
%   amplitude_threshold: Narrow band amplitude threshold (max - min)
%   max_refinement_iter: Maximum number of fine scanning iterations
%   refinement_margin: Boundary expansion ratio per iteration
%   refinement_points_factor: Points multiplication factor per iteration
%   rho: Density
%   c: Sound speed
%   verbose: Whether to display detailed information (optional, default true)
%
% Returns:
%   final_result: Final narrow band result [k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude]
%   refinement_history: Iteration history (each row is result of one iteration)
%   candidates: All candidate narrow bands
%   k_coarse: k values from coarse scan
%   real_coarse: Real part of eigenvalues from coarse scan

% Default parameters
if nargin < 18
    verbose = true;
end

%% ====== Coarse Scan ======
if verbose
    fprintf('Coarse scanning...\n');
end
[real_coarse, ~, k_coarse] = run_scan_from_seq(k_min0, k_max0, num_k_coarse, S_seq, d, rho, c);

%% ====== Slope Detection for Narrow Bands ======
if verbose
    fprintf('Detecting narrow bands with slope absolute value > %.3f and Δk/k0 < %.2f and amplitude >= %.2f...\n', ...
        slope_threshold, width_threshold, amplitude_threshold);
end

q = real_coarse(1,:);
dq = diff(q);
dk = diff(k_coarse);
slope = abs(dq ./ dk);

% Find high slope regions
high_slope_mask = slope > slope_threshold;
high_slope_idx = find(high_slope_mask);

if isempty(high_slope_idx)
    error('No region found with slope absolute value greater than %.3f', slope_threshold);
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
            candidates = [candidates; k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude, i];
        elseif verbose
            fprintf('  Excluding candidate %d: Δk/k0=%.6e (require<%.2f), amplitude=%.6e (require>=%.2f)\n', ...
                i, delta_k_norm, width_threshold, amplitude, amplitude_threshold);
        end
    end
end

if isempty(candidates)
    error('No narrow band found satisfying slope, width, and amplitude conditions simultaneously');
end

if verbose
    fprintf('Found %d narrow bands meeting conditions:\n', size(candidates,1));
    for ii = 1:size(candidates,1)
        fprintf('  Narrow band %d: k∈[%.6e, %.6e], Δk/k0=%.6e, average slope=%.6e, amplitude=%.6e\n', ...
            ii, candidates(ii,1), candidates(ii,2), candidates(ii,3), candidates(ii,4), candidates(ii,6));
    end
end

%% ====== Select Optimal Narrow Band ======
% Sort by average slope (higher slope prioritized)
[~, sort_idx] = sort(candidates(:,4), 'descend');
sorted_candidates = candidates(sort_idx, :);

if verbose
    fprintf('\nSorted by average slope:\n');
    for ii = 1:size(sorted_candidates,1)
        fprintf('  %d: k∈[%.6e, %.6e], average slope=%.6e, Δk/k0=%.6e, amplitude=%.6e\n', ...
            ii, sorted_candidates(ii,1), sorted_candidates(ii,2), ...
            sorted_candidates(ii,4), sorted_candidates(ii,3), sorted_candidates(ii,6));
    end
end

% Select narrow band with maximum average slope
best_candidate = sorted_candidates(1, :);

if verbose
    fprintf('\nSelected narrow band with maximum average slope:\n');
    fprintf('  k∈[%.8e, %.8e], Δk/k0=%.6e, average slope=%.6e, amplitude=%.6e\n', ...
        best_candidate(1), best_candidate(2), best_candidate(3), best_candidate(4), best_candidate(6));
end

%% ====== Multiple Refinement Scans ======
if verbose
    fprintf('\nStarting multiple refinement scans...\n');
end

current_candidate = best_candidate;
refinement_history = [current_candidate]; % Record result of each iteration

for iter = 1:max_refinement_iter
    if verbose
        fprintf('\n--- Refinement scan iteration %d/%d ---\n', iter, max_refinement_iter);
    end
    
    % Calculate current scanning parameters
    current_width = current_candidate(2) - current_candidate(1);
    current_margin = refinement_margin * current_width;
    current_num_points = round(num_k_fine * (refinement_points_factor^(iter-1)));
    
    % Determine scanning interval (gradually expand)
    kf_min = max(current_candidate(1) - current_margin, k_min0);
    kf_max = min(current_candidate(2) + current_margin, k_max0);
    
    if verbose
        fprintf('  Scanning interval: [%.8e, %.8e]\n', kf_min, kf_max);
        fprintf('  Number of scanning points: %d\n', current_num_points);
    end
    
    % Execute fine scan
    [real_fine, ~, k_fine] = run_scan_from_seq(kf_min, kf_max, current_num_points, S_seq, d, rho, c);
    
    % Re-detect narrow bands on fine scan data
    fine_candidates = detect_narrow_bands_by_slope(k_fine, real_fine, k0, slope_threshold, width_threshold, amplitude_threshold);
    
    if isempty(fine_candidates)
        if verbose
            fprintf('  Iteration %d: No narrow band detected meeting conditions, keeping previous result\n', iter);
        end
        break;
    end
    
    % Select narrow band overlapping with current candidate and having maximum slope
    best_fine_candidate = select_best_overlapping_candidate(fine_candidates, current_candidate);
    
    if isempty(best_fine_candidate)
        if verbose
            fprintf('  Iteration %d: No overlapping narrow band found, keeping previous result\n', iter);
        end
        break;
    end
    
    % Update current candidate
    prev_candidate = current_candidate;
    current_candidate = best_fine_candidate;
    refinement_history = [refinement_history; current_candidate];
    
    if verbose
        fprintf('  Iteration %d result: k∈[%.8e, %.8e], Δk/k0=%.6e, average slope=%.6e, amplitude=%.6e\n', ...
            iter, current_candidate(1), current_candidate(2), current_candidate(3), current_candidate(4), current_candidate(6));
    end
    
    % Check convergence conditions
    width_change = abs(current_candidate(3) - prev_candidate(3)) / prev_candidate(3);
    position_change = max(abs(current_candidate(1) - prev_candidate(1)), ...
                         abs(current_candidate(2) - prev_candidate(2))) / (prev_candidate(2) - prev_candidate(1));
    
    if verbose
        fprintf('  Width change: %.2f%%, Position change: %.2f%%\n', width_change*100, position_change*100);
    end
    
    if width_change < 0.01 && position_change < 0.02  % Change less than 1% and 2%
        if verbose
            fprintf('  Convergence condition reached, stopping iteration\n');
        end
        break;
    end
end

final_result = current_candidate;

%% ====== Final Results ======
if verbose
    fprintf('\n=== Final Results ===\n');
    fprintf('k_min = %.12e\n', final_result(1));
    fprintf('k_max = %.12e\n', final_result(2));
    fprintf('Δk/k0 = %.12e\n', final_result(3));
    fprintf('Average slope = %.12e\n', final_result(4));
    fprintf('Maximum slope = %.12e\n', final_result(5));
    fprintf('Amplitude = %.12e\n', final_result(6));
end

end

