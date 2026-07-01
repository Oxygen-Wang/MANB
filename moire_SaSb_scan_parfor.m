% moire_SaSb_scan_parfor.m
% Parallel heatmap scan of Sa,Sb with fixed N1,N2
% Using standard method find_narrow_band for narrow band detection
% Note: Sa=Sa1, Sb=Sb1
clear; clc; close all;

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
addpath(fullfile(figure_code_dir, 'functions'));  % Use local functions folder

%% ========== Parameters (Modifiable) ==========
% Fixed parameters
N1 = 12;                % Subsection length
N2 = N1 + 1;
d = lcm(N1, N2);
num_n = d * 8;
k0 = 2*pi / d;
rho = 1000; c = 1000;

dd=0.1;
% Sa,Sb scanning range
Sa_range = 1:dd:4;
Sb_range = 1:dd:4;



% Scanning parameters
k_min0 = 0;
k_max0 = 11.1*k0;
num_k_coarse = 1001;
num_k_fine = 201;

% Detection thresholds
slope_threshold = pi;       % Absolute value threshold for slope
width_threshold = 0.5;      % Δk/k0 threshold
amplitude_threshold = 1;    % Narrow band amplitude threshold

% Multiple refinement scanning parameters
max_refinement_iter =25;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

%% ========== Initialization ==========
n1 = length(Sa_range);
n2 = length(Sb_range);
total_iterations = n1 * n2;

% Pre-allocate result containers
delta_k_lin = nan(total_iterations,1);
avg_slope_lin = nan(total_iterations,1);
convergence_lin = nan(total_iterations,1);
errors = repmat({''}, total_iterations, 1);

fprintf('Starting parallel SaSb scan (N1=%d, N2=%d, total %d tasks)...\n', N1, N2, total_iterations);
fprintf('Sa range: %.2f-%.2f, Sb range: %.2f-%.2f\n', min(Sa_range), max(Sa_range), min(Sb_range), max(Sb_range));

% Start parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool('local');
end

%% ========== Parallel Main Loop ==========
parfor idx = 1:total_iterations
    try
        % Calculate i,j
        i = floor((idx-1) / n2) + 1;
        j = mod((idx-1), n2) + 1;
        Sa_val = Sa_range(i);
        Sb_val = Sb_range(j);
        
        % ====== Use standard method find_narrow_band to find narrow band ======
               [final_result, refinement_history, ~, ~, ~] = find_narrow_band(...
            N1, N2, Sa_val, Sb_val, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            false);  % verbose = false to avoid parallel output confusion
        
        % Extract results
        % final_result: [k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude]
        delta_k_lin(idx) = final_result(3);
        avg_slope_lin(idx) = final_result(4);
        convergence_lin(idx) = size(refinement_history, 1);
        
    catch ME
        % Record error message (find_narrow_band throws error when no narrow band is found)
        errors{idx} = ME.message;
        delta_k_lin(idx) = NaN;
        avg_slope_lin(idx) = NaN;
        convergence_lin(idx) = -1;
    end
end % parfor

%% ========== Restore Matrix from Linear Array ==========
delta_k_data = nan(n1, n2);
avg_slope_data = nan(n1, n2);
convergence_data = nan(n1, n2);

for idx = 1:total_iterations
    i = floor((idx-1) / n2) + 1;
    j = mod((idx-1), n2) + 1;
    delta_k_data(i,j) = delta_k_lin(idx);
    avg_slope_data(i,j) = avg_slope_lin(idx);
    convergence_data(i,j) = convergence_lin(idx);
end

%% ========== Save Results ==========
% Set output directory
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = fullfile(figure_code_dir, 'output');  % If not exists, save to figure_code/output
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
end

result_filename = fullfile(output_dir, sprintf('moire_results_SaSb_N1%d_N2%d.csv', N1, N2));
fid = fopen(result_filename, 'w');
fprintf(fid, 'Sa,Sb,Delta_k_k0,Avg_Slope,Convergence_Iterations,Error\n');
for i = 1:n1
    for j = 1:n2
        idx = (i-1)*n2 + j;
        fprintf(fid, '%.3f,%.3f,%.6e,%.6e,%d,"%s"\n', ...
            Sa_range(i), Sb_range(j), delta_k_data(i,j), avg_slope_data(i,j), convergence_data(i,j), errors{idx});
    end
end
fclose(fid);
fprintf('Results saved to: %s\n', result_filename);

%% ========== 可视化热力图 ==========
%% ===== 简化版热力图可视化 =====
figure('Position', [100, 100, 1200, 900]);

% 1. Narrow band width heatmap
subplot(2,2,1);
delta_k_plot = abs(delta_k_data);
delta_k_plot(delta_k_plot <= 0 | delta_k_plot > 0.5) = NaN;
pcolor(Sb_range, Sa_range, delta_k_plot);
shading flat; colorbar;
xlabel('Sb'); ylabel('Sa');
title('Narrow Band Width Heatmap');

% 2. Average slope heatmap
subplot(2,2,2);
pcolor(Sb_range, Sa_range, avg_slope_data);
shading flat; colorbar;
xlabel('Sb'); ylabel('Sa');
title('Average Slope Heatmap');

% 3. Convergence iteration count heatmap
subplot(2,2,3);
pcolor(Sb_range, Sa_range, convergence_data);
shading flat; colorbar;
xlabel('Sb'); ylabel('Sa');
title('Convergence Iteration Count');

% 4. Valid data points heatmap
subplot(2,2,4);
valid_mask = ~isnan(delta_k_data);
pcolor(Sb_range, Sa_range, double(valid_mask));
shading flat;
colormap(gca, [1 1 1; 0 0.5 0]); % White-Green
colorbar('Ticks', [0 1], 'TickLabels', {'No narrow band', 'Has narrow band'});
xlabel('Sb'); ylabel('Sa');
title('Valid Narrow Band Detection');

% Save image
saveas(gcf, fullfile(output_dir, sprintf('moire_heatmap_N1%d_N2%d_SaSb_step%.1f.png', N1, N2, dd)));

%% ========== Statistics ==========
fprintf('\n=== Statistics ===\n');
fprintf('Total scanning points: %d\n', total_iterations);
fprintf('Points with detected narrow bands: %d (%.1f%%)\n', sum(~isnan(delta_k_data(:))), 100*sum(~isnan(delta_k_data(:)))/total_iterations);
fprintf('Narrow band width range: %.3e - %.3e\n', min(delta_k_data(:), [], 'omitnan'), max(delta_k_data(:), [], 'omitnan'));
fprintf('Average slope range: %.3e - %.3e\n', min(avg_slope_data(:), [], 'omitnan'), max(avg_slope_data(:), [], 'omitnan'));

%% ====== Find Optimal Position ======
% Remove NaN
valid_mask = ~isnan(delta_k_data);

% If you want to find minimum narrow band width
[min_val, idx] = min(delta_k_data(valid_mask));

% Convert idx back to 2D indices
[row, col] = ind2sub(size(delta_k_data), find(valid_mask));
best_row = row(idx);
best_col = col(idx);

% Output optimal Sa, Sb
Sa_best = Sa_range(best_row);
Sb_best = Sb_range(best_col);

fprintf('\n=== Optimal Solution (Based on Minimum Narrow Band Width) ===\n');
fprintf('Minimum narrow band width = %.3e\n', min_val);
fprintf('Best Sa = %.4f (Sa1 = Sa = %.4f)\n', Sa_best, Sa_best);
fprintf('Best Sb = %.4f (Sb1 = Sb = %.4f)\n', Sb_best, Sb_best);







% Note: All helper functions (run_scan, gener_n, acoustic_TM, detect_narrow_bands_by_slope, 
% select_best_overlapping_candidate, contiguous_groups, etc.) are in the functions directory,
% called by find_narrow_band function, no need to redefine in this file.