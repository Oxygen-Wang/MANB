% Example script using find_narrow_band function to find narrow bands
% 
% This script demonstrates how to directly input parameters to solve for system narrow bands

clear; clc;

% Add functions folder to path (if needed)
if exist('functions', 'dir')
    addpath('functions');
end

%% ========== Parameter Settings ==========

% System parameters
N1 = 12;
N2 = 13;
d1 = 4;
d2 = 4;

% Scanning parameters
k_min0 = 0;
k_max0 = 0.4541;
num_k_coarse = 1001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter = 3;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Physical parameters
rho = 1000;  % Density
c = 1000;    % Sound speed

%% ========== Call Function to Find Narrow Band ==========

fprintf('Starting to find narrow band...\n\n');

[final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
    find_narrow_band(...
    N1, N2, d1, d2, ...
    k_min0, k_max0, ...
    num_k_coarse, num_k_fine, ...
    slope_threshold, width_threshold, amplitude_threshold, ...
    max_refinement_iter, refinement_margin, refinement_points_factor, ...
    rho, c, ...
    true);  % verbose = true to show detailed information

%% ========== Result Output ==========

fprintf('\n\n=== Final Narrow Band Results ===\n');
fprintf('k_min = %.12e\n', final_result(1));
fprintf('k_max = %.12e\n', final_result(2));
fprintf('Δk/k0 = %.12e\n', final_result(3));
fprintf('Average slope = %.12e\n', final_result(4));
fprintf('Maximum slope = %.12e\n', final_result(5));
fprintf('Amplitude = %.12e\n', final_result(6));

%% ========== Plotting (consistent with original) ==========

% Calculate slope data (for plotting)
q = real_coarse(1,:);
dq = diff(q);
dk = diff(k_coarse);
slope = abs(dq ./ dk);

% Sort candidates by average slope (for plotting)
[~, sort_idx] = sort(candidates(:,4), 'descend');
sorted_candidates = candidates(sort_idx, :);

% Calculate basic parameters (for final precision scan)
d = lcm(N1, N2);
num_n = d * 8;
k0 = 2*pi / d;
Sa = d1;
Sb = d2;
Sa1 = d1;
Sb1 = d2;
n1_t = gener_n(num_n, N2, Sa, Sb);
n2_t = gener_n(num_n, N1, Sa1, Sb1);
n_pro = n1_t + n2_t;

% Execute final precision scan for plotting
kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
[real_final, ~, k_final] = run_scan(kf_final_min, kf_final_max, num_k_fine*2, n_pro, d, num_n, rho, c);

% Calculate slope of final narrow band
q_final = real_final(1,:);
dq_final = diff(q_final);
dk_final = diff(k_final);
slope_final = abs(dq_final ./ dk_final);

figure('Position', [100, 100, 1400, 1000]);

% Main plot: Band structure and narrow band evolution
subplot(2,3,1);
plot(k_coarse, real_coarse(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(k_coarse, real_coarse(2,:), 'r-', 'LineWidth', 1.5);
% Mark all candidate narrow bands
best_candidate_idx = sort_idx(1);  % Index of optimal narrow band in original candidates
for i = 1:size(candidates,1)
    color = [0.9, 0.9, 0.9];
    if i == best_candidate_idx
        color = [0.8, 1, 0.8];
    end
    fill([candidates(i,1), candidates(i,2), candidates(i,2), candidates(i,1)], ...
         [min(real_coarse(:)), min(real_coarse(:)), max(real_coarse(:)), max(real_coarse(:))], ...
         color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
% Mark iteration history
colors = lines(size(refinement_history,1));
for i = 1:size(refinement_history,1)
    plot([refinement_history(i,1), refinement_history(i,1)], [min(real_coarse(:)), max(real_coarse(:))], ...
         '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot([refinement_history(i,2), refinement_history(i,2)], [min(real_coarse(:)), max(real_coarse(:))], ...
         '--', 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlabel('Wavenumber k'); ylabel('Bloch phase q');
title('Band Structure and Narrow Band Evolution');
legend('Band 1', 'Band 2', 'Candidate narrow bands', 'Optimal narrow band', 'Iteration boundaries', 'Location', 'best');
grid on;

% Subplot 2: Slope analysis
subplot(2,3,2);
plot(k_coarse(1:end-1), slope, 'k-', 'LineWidth', 1.5); hold on;
plot([k_min0, k_max0], [slope_threshold, slope_threshold], 'r--', 'LineWidth', 1.5);
% Mark slopes within candidate narrow bands
for i = 1:size(candidates,1)
    idx_start = find(k_coarse >= candidates(i,1), 1);
    idx_end = find(k_coarse <= candidates(i,2), 1, 'last') - 1;
    if ~isempty(idx_start) && ~isempty(idx_end) && idx_end >= idx_start
        color = [0.7, 0.7, 0.7];
        if i == best_candidate_idx
            color = [0, 0.8, 0];
        end
        plot(k_coarse(idx_start:idx_end), slope(idx_start:idx_end), '-', 'Color', color, 'LineWidth', 2);
    end
end
xlabel('Wavenumber k'); ylabel('Slope |dq/dk|');
title(sprintf('Band Slope (threshold=%.3f)', slope_threshold));
legend('Slope', 'Threshold', 'Candidate narrow bands', 'Optimal narrow band', 'Location', 'best');
grid on;

% Subplot 3: Iteration history comparison
subplot(2,3,3);
if size(refinement_history,1) > 1
    plot(0:size(refinement_history,1)-1, refinement_history(:,3), 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Iteration number'); ylabel('Δk/k0');
    title('Narrow Band Width Evolution');
    grid on;
end

% Subplot 4: Final precision scan results
subplot(2,3,4);
plot(k_final, real_final(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(k_final, real_final(2,:), 'r-', 'LineWidth', 1.5);
% Mark final narrow band
ymin = min(real_final(:)); ymax = max(real_final(:));
fill([final_result(1), final_result(2), final_result(2), final_result(1)], ...
     [ymin, ymin, ymax, ymax], [0.8, 1, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'LineWidth', 2);
xlabel('Wavenumber k'); ylabel('Bloch phase q');
title('Final Precision Scan Results');
legend('Band 1', 'Band 2', 'Final narrow band', 'Location', 'best');
grid on;

% Subplot 5: Final narrow band slope analysis
subplot(2,3,5);
plot(k_final(1:end-1), slope_final, 'k-', 'LineWidth', 1.5); hold on;
% Mark slopes within final narrow band
idx_start_final = find(k_final >= final_result(1), 1);
idx_end_final = find(k_final <= final_result(2), 1, 'last') - 1;
if ~isempty(idx_start_final) && ~isempty(idx_end_final) && idx_end_final >= idx_start_final
    plot(k_final(idx_start_final:idx_end_final), slope_final(idx_start_final:idx_end_final), ...
         'r-', 'LineWidth', 2);
    plot([final_result(1), final_result(1)], [0, max(slope_final)], 'g--', 'LineWidth', 1.5);
    plot([final_result(2), final_result(2)], [0, max(slope_final)], 'g--', 'LineWidth', 1.5);
end
xlabel('Wavenumber k'); ylabel('Slope |dq/dk|');
title('Final Narrow Band Slope Analysis');
legend('Overall slope', 'Slope within narrow band', 'Narrow band boundaries', 'Location', 'best');
grid on;

% Subplot 6: Candidate narrow band comparison
subplot(2,3,6);
if size(candidates,1) > 1
    scatter(candidates(:,3), candidates(:,4), 80, 'b', 'filled'); hold on;
    scatter(final_result(3), final_result(4), 120, 'r', 'filled', '^');
    
    xlabel('Δk/k0'); ylabel('Average slope');
    title('Candidate Narrow Band Comparison');
    legend('Candidates', 'Final', 'Location', 'best');
    grid on;
end


