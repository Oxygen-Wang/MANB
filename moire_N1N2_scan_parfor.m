%% moire_N1N2_scan_parfor.m
% Parallelized N1,N2 scan using parfor - using standard method find_narrow_band
clear; clc; close all;

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
addpath(fullfile(figure_code_dir, 'functions'));  % Use local functions folder

tic
%% ========== Modifiable Parameters ==========
% N1,N2 scanning range
N1_range = 1:1:20; % example: 1..20
N2_range = 1:1:20; % example: 1..20

% Fixed parameters
d1 = 1.0000
d2= 2
rho = 1000; 
c = 1000;

% Multiple of k_max0
k_max_multiple = 11;  % k_max0 = 11 * k0

% Scanning parameters
k_min0 = 0; 
num_k_coarse = 2001;
num_k_fine =500;

% Detection thresholds
slope_threshold = pi / 0.5;  % Absolute value threshold for slope
width_threshold = 1;         % Δk/k0 threshold
amplitude_threshold = 1;     % New: narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =25;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

%% ========== Initialization ==========
n1 = length(N1_range);
n2 = length(N2_range);
total_iterations = n1 * n2;

delta_k_lin = nan(total_iterations,1);
avg_slope_lin = nan(total_iterations,1);
convergence_lin = nan(total_iterations,1);
errors = repmat({''}, total_iterations, 1);

fprintf('Starting parallel N1N2 scan (total %d tasks)...\n', total_iterations);
fprintf('Detection thresholds: slope>%.3f, width<%.2f, amplitude>=%.2f\n', ...
    slope_threshold, width_threshold, amplitude_threshold);
fprintf('k_max0 = %d * k0\n', k_max_multiple);

%% ========== Parallel Main Loop ==========
parfor idx = 1:total_iterations
    try
        % Calculate current N1,N2 indices
        i = floor((idx-1) / n2) + 1;
        j = mod((idx-1), n2) + 1;
        N1 = N1_range(i);
        N2 = N2_range(j);

        % Calculate period d of current structure (for calculating k0)
        d = lcm(N1, N2);
        k0 = 2 * pi / d;  % Each N1,N2 combination has corresponding k0
        
        % Dynamically calculate k_max0
        k_max0 = k_max_multiple * k0;

        % ====== Use standard method find_narrow_band to find narrow band ======
        % Call find_narrow_band function (verbose=false to avoid parallel output confusion)
        [final_result, refinement_history, ~, ~, ~] = find_narrow_band(...
            N1, N2, d1, d2, ...
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
end

%% ========== Restore Matrix and Save ==========
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

%% === Modified Section: Safe File Writing ===
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

result_filename = fullfile(output_dir, 'moire_results_N1N2_parfor_with_amplitude.csv');
[fid, msg] = fopen(result_filename, 'w');
if fid == -1
    error('Cannot open output file: %s\nError message: %s', result_filename, msg);
end

fprintf(fid, 'N1,N2,Delta_k_k0,Avg_Slope,Convergence_Iterations,Error\n');
for i = 1:n1
    for j = 1:n2
        idx = (i-1)*n2 + j;
        fprintf(fid, '%d,%d,%.6e,%.6e,%d,"%s"\n', ...
            N1_range(i), N2_range(j), delta_k_data(i,j), ...
            avg_slope_data(i,j), convergence_data(i,j), errors{idx});
    end
end
fclose(fid);
fprintf('Results saved to: %s\n', result_filename);

%% ========== Visualize Heatmaps ==========
Delta_k_plot = abs(delta_k_data);
Delta_k_plot(Delta_k_plot<=0) = NaN;
Delta_k_plot(Delta_k_plot>0.5) = NaN;

figure('Position',[100,100,1000,800]);
h = pcolor(N2_range, N1_range, Delta_k_plot);
set(h,'EdgeColor','none');
set(gca,'ColorScale','log');
colormap(jet); colorbar;
xlabel('N2'); ylabel('N1');
title(sprintf('Narrow Band Width Heatmap (\\Delta k / k_0) [log scale] - k_max0 = %d*k0', k_max_multiple));
set(gca,'XTick', N2_range);
set(gca,'YTick', N1_range);
saveas(gcf, fullfile(output_dir,'moire_heatmap_N1N2_log_parfor_with_amplitude.png'));

figure('Position',[100,100,1000,800]);
h2 = pcolor(N2_range, N1_range, avg_slope_data);
set(h2,'EdgeColor','none');
colormap(jet); colorbar;
xlabel('N2'); ylabel('N1');
title(sprintf('Average Slope Heatmap - k_max0 = %d*k0', k_max_multiple));
set(gca,'XTick', N2_range);
set(gca,'YTick', N1_range);
saveas(gcf, fullfile(output_dir,'moire_heatmap_N1N2_slope_parfor_with_amplitude.png'));
toc

%% ========== Print Table ==========
fprintf('\n=== Results Summary (Δk/k0) ===\n');
fprintf('N1\\N2');
for j = N2_range, fprintf('\t%d', j); end
fprintf('\n');
for i = 1:n1
    fprintf('%d', N1_range(i));
    for j = 1:n2
        if ~isnan(delta_k_data(i,j))
            fprintf('\t%.3e', delta_k_data(i,j));
        else
            fprintf('\t---');
        end
    end
    fprintf('\n');
end

% Note: All helper functions (run_scan, gener_n, acoustic_TM, detect_narrow_bands_by_slope, 
% select_best_overlapping_candidate, contiguous_groups, etc.) are in the functions directory,
% called by find_narrow_band function, no need to redefine in this file.
