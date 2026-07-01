clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'functions'));
%% ========== Parameter Settings ==========
% N1 value list
N1_list = [4, 12,16];
d1 = 1;
d2 = 4;
% Scanning parameters
k_min0 = 0;
num_k_coarse = 5001;
num_k_fine =200;
% Detection thresholds
slope_threshold = pi;
width_threshold = 0.5;
amplitude_threshold = 1;
% Multiple refinement scanning parameters
max_refinement_iter = 25;
refinement_margin = 0.3;
refinement_points_factor = 1.5;
% Physical parameters
rho = 1000;  % Density
c = 1000;    % Sound speed
% Sensitivity calculation parameters
n_d =1;               % Number of periods (for transmission calculation)
num_k_sens =5001;       % Number of frequency scanning points
%% ========== Store Results ==========
results = cell(length(N1_list), 1);
sensitivity_max_list = zeros(length(N1_list), 1);
sensitivity_peak_delta_k_list = zeros(length(N1_list), 1);
k_star_list = zeros(length(N1_list), 1);
k_end_list = zeros(length(N1_list), 1);
%% ========== Calculate for Each N1 Value ==========
for idx = 1:length(N1_list)
    N1 = N1_list(idx);
    N2 = N1 + 1;  % N2 = N1 + 1
    fprintf('\n\n========================================\n');
    fprintf('Calculating case N1 = %d, N2 = %d\n', N1, N2);
    fprintf('========================================\n\n');
    % Calculate basic parameters
    d = lcm(N1, N2);
    k0 = 2*pi/d;
    k_max0 = 11.1*k0;
    % ========== Find Narrow Band ==========
    fprintf('Starting to find narrow band...\n\n');
    [final_result, ~, ~, ~, ~] = ...
        find_narrow_band(...
        N1, N2, d1, d2, ...
        k_min0, k_max0, ...
        num_k_coarse, num_k_fine, ...
        slope_threshold, width_threshold, amplitude_threshold, ...
        max_refinement_iter, refinement_margin, refinement_points_factor, ...
        rho, c, ...
        false);  % verbose = false to reduce output

    if isempty(final_result)
        fprintf('Warning: No narrow band found for N1 = %d, skipping\n', N1);
        continue;
    end

    % Calculate center frequency
    k_star_list(idx) =final_result(1);
    k_end_list(idx)=final_result(2);  
    % ========== Prepare Structure Parameters ==========
    num_n = d * 8;
    Sa = d1; Sb = d2;
    Sa1 = d1; Sb1 = d2;
    n1_t = gener_n(num_n, N2, Sa, Sb);
    n2_t = gener_n(num_n, N1, Sa1, Sb1);
    n_pro = n1_t + n2_t;
    k0 = 2*pi/d;
    % ========== Calculate Sensitivity ==========
    fprintf('\n=== Starting Sensitivity Calculation ===\n');
    if N1==N1_list(1)
        [sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity(...
            k_star_list(1), n_pro, d, num_n, rho, c, n_d, k_end_list(1), num_k_sens, ...
            [], false);  % Don't save separate figure
    elseif N1==N1_list(2)
        [sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity(...
            k_star_list(2), n_pro, d, num_n, rho, c, n_d,  k_end_list(2), num_k_sens, ...
            [], false);  % Don't save separate figure
    elseif N1==N1_list(3)
        [sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity(...
            k_star_list(3), n_pro, d, num_n, rho, c, n_d, k_end_list(3), num_k_sens, ...
            [], false);  % Don't save separate figure
    end
    % Store results
    results{idx} = struct();
    results{idx}.N1 = N1;
    results{idx}.N2 = N2;  % N2 = N1 + 1
    results{idx}.d = d;
    results{idx}.k0 = k0;
    results{idx}.sensitivity_result = sensitivity_result;
    results{idx}.delta_k_values = delta_k_values;
    results{idx}.intensity2 = intensity2;
    
    sensitivity_max_list(idx) = sensitivity_result.sensitivity_max;
    sensitivity_peak_delta_k_list(idx) = sensitivity_result.sensitivity_peak_delta_k;
    
    fprintf('Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
    fprintf('Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
    
    % Calculate sensitivity order of magnitude (power of 10)
    sensitivity_order = floor(log10(sensitivity_result.sensitivity_max));
    fprintf('Sensitivity order of magnitude: 10^%d\n', sensitivity_order);
end

%% ========== Plot Three Subplots ==========
figure('Position', [100, 100, 1400, 500]);

for idx = 1:length(N1_list)
    if isempty(results{idx})
        continue;
    end
    
    subplot(1, 3, idx);
    
    % =========================
    % Left axis: log intensity
    % =========================
    yyaxis left
    log_intensity = log(results{idx}.intensity2 / min(results{idx}.intensity2));
    plot(results{idx}.delta_k_values, log_intensity, 'b-', 'LineWidth', 2.5);
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    ylim([min(log_intensity), max(log_intensity)]);  % ⭐关键
    
    % =========================
    % Right axis: sensitivity
    % =========================
    yyaxis right
    sens = results{idx}.sensitivity_result.sensitivity;
    plot(results{idx}.delta_k_values, sens, 'r-', 'LineWidth', 2.5);
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    ylim([min(sens), max(sens)]);  % ⭐关键
    
    % =========================
    % X axis
    % =========================
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    xlim([min(results{idx}.delta_k_values), max(results{idx}.delta_k_values)]); % ⭐关键
    
    % =========================
    % Title & style
    % =========================
    title(sprintf('$N_1 = %d$', results{idx}.N1), ...
          'Interpreter', 'latex', 'FontSize', 18);
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

sgtitle('Moiré Structure Sensitivity Comparison', ...
        'FontSize', 22, 'FontWeight', 'bold');

%% ========== Save Figures ==========
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end

% Save EPS format
print('-depsc', fullfile(output_dir, 'Fig_Moire_Sensitivity_N1_Comparison.eps'));
% Save PNG format (high resolution)
print('-dpng', '-r300', fullfile(output_dir, 'Fig_Moire_Sensitivity_N1_Comparison.png'));

fprintf('\n\nImages saved to: %s\n', fullfile(output_dir, 'Fig_Moire_Sensitivity_N1_Comparison.eps'));
fprintf('Images saved to: %s\n', fullfile(output_dir, 'Fig_Moire_Sensitivity_N1_Comparison.png'));

%% ========== Export Sensitivity Data to CSV ==========
fprintf('\n\n=== Starting to Export Sensitivity Data to CSV ===\n');

% Create figure directory
script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

% Export sensitivity data for each N1 value
for idx = 1:length(N1_list)
    if isempty(results{idx})
        fprintf('Skipping N1 = %d: No data\n', N1_list(idx));
        continue;
    end
    
    delta_k_values = results{idx}.delta_k_values(:);
    intensity2 = results{idx}.intensity2(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = results{idx}.sensitivity_result.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    if length(sens) < length(delta_k_values)
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    

% 从最大灵敏度位置开始保存
[~, peak_idx] = max(abs(sens));

delta_k_values = delta_k_values(peak_idx:end);
intensity2     = intensity2(peak_idx:end);
sens           = sens(peak_idx:end);
log_intensity  = log_intensity(peak_idx:end);

    
    % Take maximum length and pad with NaN
    Ls = max([length(delta_k_values), length(intensity2), length(sens), length(log_intensity)]);
    dk_pad = nan(Ls,1); 
    dk_pad(1:length(delta_k_values)) = delta_k_values;
    
    int_pad = nan(Ls,1); 
    int_pad(1:length(intensity2)) = intensity2;
    
    sens_pad = nan(Ls,1); 
    sens_pad(1:length(sens)) = sens;
    
    log_int_pad = nan(Ls,1);
    log_int_pad(1:length(log_intensity)) = log_intensity;
    
    % Save CSV file
    sens_file = fullfile(figure_dir, sprintf('Sensitivity_Moire_N1_%d.csv', results{idx}.N1));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end

fprintf('\n=== Sensitivity Data Export Completed ===\n');

%% ========== Output Results Summary ==========
fprintf('\n\n========================================\n');
fprintf('Sensitivity Calculation Results Summary\n');
fprintf('========================================\n\n');

for idx = 1:length(N1_list)
    if isempty(results{idx})
        continue;
    end
    
    N1 = results{idx}.N1;
    N2 = results{idx}.N2;
    sensitivity_max = sensitivity_max_list(idx);
    sensitivity_order = floor(log10(sensitivity_max));
    
    fprintf('N1 = %2d, N2 = %2d:  Maximum sensitivity = %.6e  (10^%d order of magnitude)\n', ...
            N1, N2, sensitivity_max, sensitivity_order);
    fprintf('                     Peak position: δk = %.6e\n', sensitivity_peak_delta_k_list(idx));
end

fprintf('========================================\n');
fprintf('Sensitivity increasing trend with N1:\n');
for idx = 1:length(N1_list)
    if ~isempty(results{idx})
        sensitivity_order = floor(log10(sensitivity_max_list(idx)));
        fprintf('N1 = %d, N2 = %d: 10^%d\n', results{idx}.N1, results{idx}.N2, sensitivity_order);
    end
end
fprintf('========================================\n');

fprintf('\n=== All Calculations Completed ===\n');

