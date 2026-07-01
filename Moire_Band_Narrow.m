% Moiré band structure and narrow band marking
% Plot band structure and narrow band marking, including main plot and top-left subplot

clear; clc; close all;

% Add functions directory to path
% script_dir = fileparts(mfilename('fullpath'));
% figure_code_dir = fileparts(script_dir);  % figure_code folder
% addpath(fullfile(figure_code_dir, 'functions'));  % Use local functions folder

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'functions'));
%% ========== Parameter Settings ==========

% System parameters
N1 =16;
N2 =17;
% N1 = 12;
% N2 =13;
% N1 = 16;
% N2 =17;
d1= 1;
d2= 4;
d=lcm(N1,N2);
k0=2*pi/d;
% Scanning parameters
k_min0 = 0;
%k_max0 = 0.4541;
k_max0 = 11.1*k0;

num_k_coarse = 1001;
num_k_fine = 200;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter = 25;     % Maximum refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Physical parameters
rho = 1000;  % Density
c = 1000;    % Sound speed

% ========== Call function to find narrow band ==========

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

% ========== Result Output ==========

fprintf('\n\n=== Final Narrow Band Results ===\n');
fprintf('k_min = %.12e\n', final_result(1));
fprintf('k_max = %.12e\n', final_result(2));
fprintf('Δk/k0 = %.12e\n', final_result(3));
fprintf('k_min/k0 = %.12e\n', final_result(1)/k0);
fprintf('k_max/k0 = %.12e\n', final_result(2)/k0);
fprintf('(k_max+k_min)/(2*k0) = %.12e\n', (final_result(2)+final_result(1))/(2*k0));
fprintf('Average slope = %.12e\n', final_result(4));
fprintf('Maximum slope = %.12e\n', final_result(5));
fprintf('Amplitude = %.12e\n', final_result(6));

% ========== Calculate basic parameters (for final precision scan) ==========

d = lcm(N1, N2);
num_n = d * 8;
k0 = 2*pi / d;
Sa = d1; Sb = d2;
Sa1 = d1; Sb1 = d2;
n1_t = gener_n(num_n, N2, Sa, Sb);
n2_t = gener_n(num_n, N1, Sa1, Sb1);
n_pro = n1_t + n2_t;

% Execute final precision scan for plotting
kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
[real_final, ~, k_final] = run_scan(kf_final_min, kf_final_max, num_k_fine*2, n_pro, d, num_n, rho, c);

%% ========== Plotting ==========

figure('Position', [100, 100, 1000, 700]);

% Main plot: Band structure and narrow band evolution (Figure 1)
axes('Position', [0.15, 0.1, 0.8, 0.85]);  % Main plot position (leave more space on left for y-axis label)
plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25); hold on;
plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25);

% Sort candidates by average slope (for plotting)
[~, sort_idx] = sort(candidates(:,4), 'descend');
best_candidate_idx = sort_idx(1);  % Index of optimal narrow band in original candidates

% Mark all candidate narrow bands
y_data_min = min(real_coarse(:));
y_data_max = max(real_coarse(:));
for i = 1:size(candidates,1)
    color = [0.9, 0.9, 0.9];
    if i == best_candidate_idx
        color = [0.8, 1, 0.8];
    end
    fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
         [y_data_min, y_data_min, y_data_max, y_data_max], ...
         color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Mark iteration history
colors = lines(size(refinement_history,1));
for i = 1:size(refinement_history,1)
    plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
         '--', 'Color', colors(i,:), 'LineWidth', 2.25);
    plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
         '--', 'Color', colors(i,:), 'LineWidth', 2.25);
end

% Frame final narrow band with green box (in main plot, expand 0.25*k0 on both sides)
% Set y coordinates according to actual data range
ymin = min(real_coarse(:)); ymax = max(real_coarse(:));
rect_left = final_result(1)/k0 - 0.25;  % Expand left boundary by 0.25*k0
rect_width = (final_result(2)-final_result(1))/k0 + 0.5;  % Increase width by 0.5*k0
rectangle('Position', [rect_left, ymin, rect_width, ymax-ymin], ...
          'EdgeColor', 'g', 'LineWidth', 3, 'LineStyle', '-');

% Set axis range, remove blank space
xlim([min(k_coarse/k0), max(k_coarse/k0)]);
ylim([-pi, pi]);

% Set y-axis ticks as multiples of pi
yticks([-pi, 0, pi]);
yticklabels({'$-\pi$', '$0$', '$\pi$'});

% Set x-axis ticks
xticks([0, 5, 10]);
xticklabels({'$0$', '$5$', '$10$'});

xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 32);
ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 32);
% Adjust label positions to ensure text is within the plot
set(get(gca, 'XLabel'), 'Units', 'normalized', 'Position', [0.25, 0, 0]);  % xlabel close to x-axis but not outside plot
set(get(gca, 'YLabel'), 'Units', 'normalized', 'Position', [-0.05, 0.5, 0]);  % ylabel close to y-axis but not outside plot

grid on;
set(gca, 'FontSize', 32);
set(gca, 'TickLabelInterpreter', 'latex');

% Subplot: Final precision scan results (Figure 4, placed at top-left corner)
axes('Position', [0.15, 0.63, 0.35, 0.32]);  % Top-left subplot position (ensure within main plot range)
plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
% Set axis range, remove blank regions
xlim([min(k_final/k0), max(k_final/k0)]);
ylim([min(real_final(:)), max(real_final(:))]);
set(gca, 'XTickLabel', []);  % Remove x-axis labels
set(gca, 'YTickLabel', []);  % Remove y-axis labels
grid on;
set(gca, 'FontSize', 11);

%% ========== Save as EPS and PNG ==========
% Save to output directory
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
% Save EPS format
print('-depsc', fullfile(output_dir, 'Fig_Moire_Band_Narrow.eps'));
% Save PNG format (high resolution)
print('-dpng', '-r300', fullfile(output_dir, 'Fig_Moire_Band_Narrow.png'));

fprintf('\nImage saved to: %s\n', fullfile(output_dir, 'Fig_Moire_Band_Narrow.eps'));
fprintf('Image saved to: %s\n', fullfile(output_dir, 'Fig_Moire_Band_Narrow.png'));

%% ========== Sensitivity Calculation Parameters ==========
n_d =10;               % Number of periods (for transmission calculation)
num_k_sens =101;        % Number of frequency scanning points
%% ========== Sensitivity Calculation ==========
fprintf('\n=== Starting sensitivity calculation ===\n');
% Calculate center frequency
k_center = final_result(1);
delta_k = final_result(2);
fprintf('Narrow band center frequency k_center = %.12e (k_center/k0 = %.6f)\n', k_center, k_center/k0);



%%

% Call sensitivity calculation function
[sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity(...
    k_center, n_pro, d, num_n, rho, c, n_d, delta_k, num_k_sens, output_dir, true);

% Extract results
sensitivity_max = sensitivity_result.sensitivity_max;
sensitivity_peak_delta_k = sensitivity_result.sensitivity_peak_delta_k;

%% ========== Export data uniformly to current_folder/figure ==========

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to export all CSV data ===\n');

% ----------------------------
% 1. Export band data (coarse scan)
% ----------------------------
% Handle real_coarse shape: support 2xN or Nx2
try
    if size(real_coarse,1) == 2
        band1 = real_coarse(1,:)';
        band2 = real_coarse(2,:)';
    elseif size(real_coarse,2) == 2
        band1 = real_coarse(:,1);
        band2 = real_coarse(:,2);
    else
        % Try transpose and retry
        real_coarse_transposed = real_coarse';
        if size(real_coarse_transposed,1) == 2
            band1 = real_coarse_transposed(1,:)';
            band2 = real_coarse_transposed(2,:)';
        else
            error('Cannot recognize real_coarse dimensions (neither 2xN nor Nx2)');
        end
    end
catch ME
    warning(ME.identifier, 'Failed to parse real_coarse: %s. Skipping coarse scan band export.', ME.message);
    band1 = [];
    band2 = [];
end

if ~isempty(band1) && ~isempty(band2)
    % Align (pad with NaN to maximum length)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, 'Band_Moire_Coarse.csv');
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);
end

% ----------------------------
% 2. Export precision scan band data
% ----------------------------
if ~isempty(k_final) && ~isempty(real_final)
    % Handle real_final shape: support 2xN or Nx2
    try
        if size(real_final,1) == 2
            fband1 = real_final(1,:)';
            fband2 = real_final(2,:)';
        elseif size(real_final,2) == 2
            fband1 = real_final(:,1);
            fband2 = real_final(:,2);
        else
            real_final_transposed = real_final';
            if size(real_final_transposed,1) == 2
                fband1 = real_final_transposed(1,:)';
                fband2 = real_final_transposed(2,:)';
            else
                error('Cannot recognize real_final dimensions (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning(ME.identifier, 'Failed to parse real_final: %s. Skipping precision band export.', ME.message);
        fband1 = [];
        fband2 = [];
    end

    if ~isempty(fband1) && ~isempty(fband2)
        % Align (pad with NaN to maximum length)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, 'Band_Moire_Final.csv');
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported precision scan bands: %s\n', final_file);
    end
else
    fprintf('No precision scan data, skipping Final export\n');
end

% ----------------------------
% 3. Export sensitivity data
% ----------------------------
if ~isempty(sensitivity_result) && ~isempty(delta_k_values) && ~isempty(intensity2)
    delta_k_vals = delta_k_values(:);
    intensity_vals = intensity2(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_vals)
        delta_k_vals = delta_k_vals';
    end
    if isrow(intensity_vals)
        intensity_vals = intensity_vals';
    end
    
    % Get sensitivity data (save absolute value)
    sens = abs(sensitivity_result.sensitivity);
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is the same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_vals)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_vals) - length(sens), 1)];
    elseif length(sens) > length(delta_k_vals)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_vals));
    end
    
    % Calculate log_intensity, handle NaN values
    valid_idx = ~isnan(intensity_vals) & intensity_vals > 0;
    log_intensity = zeros(size(intensity_vals)) * NaN; % Initialize with NaN
    if any(valid_idx)
        min_intensity = min(intensity_vals(valid_idx));
        if min_intensity > 0
            log_intensity(valid_idx) = log(intensity_vals(valid_idx) / min_intensity);
        end
    end
    
    % Take maximum length and pad with NaN
    Ls = max([length(delta_k_vals), length(intensity_vals), length(sens), length(log_intensity)]);
    dk_pad = nan(Ls,1); 
    dk_pad(1:length(delta_k_vals)) = delta_k_vals;
    
    int_pad = nan(Ls,1); 
    int_pad(1:length(intensity_vals)) = intensity_vals;
    
    sens_pad = nan(Ls,1); 
    sens_pad(1:length(sens)) = sens;
    
    log_int_pad = nan(Ls,1);
    log_int_pad(1:length(log_intensity)) = log_intensity;

    sens_file = fullfile(figure_dir, 'Sensitivity_Moire.csv');
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
else
    fprintf('No sensitivity data or calculation failed, skipping sensitivity export\n');
end

fprintf('\n=== All CSV data export completed ===\n');

%% ========== Save parameters to file ==========
param_file = fullfile(output_dir, 'Moire_Narrow_Band_Parameters.txt');
fid = fopen(param_file, 'w');
if fid ~= -1
    fprintf(fid, '===== Moiré Narrow Band Detection Parameters =====\n\n');
    fprintf(fid, 'System Parameters:\n');
    fprintf(fid, '  N1 = %d\n', N1);
    fprintf(fid, '  N2 = %d\n', N2);
    fprintf(fid, '  d1 = %d\n', d1);
    fprintf(fid, '  d2 = %d\n', d2);
    fprintf(fid, '  d = lcm(N1, N2) = %d\n', d);
    fprintf(fid, '  k0 = 2π/d = %.12e\n', k0);
    fprintf(fid, '  ρ = %.0f kg/m³\n', rho);
    fprintf(fid, '  c = %.0f m/s\n\n', c);
    
    fprintf(fid, 'Scanning Parameters:\n');
    fprintf(fid, '  k_min0 = %.12e\n', k_min0);
    fprintf(fid, '  k_max0 = %.12e (%.2f k0)\n', k_max0, k_max0/k0);
    fprintf(fid, '  num_k_coarse = %d\n', num_k_coarse);
    fprintf(fid, '  num_k_fine = %d\n\n', num_k_fine);
    
    fprintf(fid, 'Detection Thresholds:\n');
    fprintf(fid, '  slope_threshold = %.3f\n', slope_threshold);
    fprintf(fid, '  width_threshold = %.2f\n', width_threshold);
    fprintf(fid, '  amplitude_threshold = %.2f\n\n', amplitude_threshold);
    
    fprintf(fid, 'Narrow Band Detection Results:\n');
    fprintf(fid, '  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
    fprintf(fid, '  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
    fprintf(fid, '  k_center = %.12e (%.6f k0)\n', k_center, k_center/k0);
    fprintf(fid, '  Δk/k0 = %.12e\n', final_result(3));
    fprintf(fid, '  Average slope = %.12e\n', final_result(4));
    fprintf(fid, '  Maximum slope = %.12e\n', final_result(5));
    fprintf(fid, '  Amplitude = %.12e\n\n', final_result(6));
    
    fprintf(fid, 'Sensitivity Calculation Parameters:\n');
    fprintf(fid, '  n_d (number of periods) = %d\n', n_d);
    fprintf(fid, '  delta_k (scanning range) = %.6e\n', delta_k);
    fprintf(fid, '  Total length L = n_d × d = %d × %d = %d\n\n', n_d, d, n_d*d);
    
    fprintf(fid, 'Sensitivity Calculation Results:\n');
    fprintf(fid, '  Maximum sensitivity = %.6e\n', sensitivity_max);
    fprintf(fid, '  Peak position δk = %.6e\n', sensitivity_peak_delta_k);
    
    fclose(fid);
    fprintf('\nParameters saved to: %s\n', param_file);
else
    fprintf('Warning: Cannot create parameter file\n');
end

fprintf('\n=== All calculations completed ===\n');





k_min= final_result(1);
k_max = final_result(2);

