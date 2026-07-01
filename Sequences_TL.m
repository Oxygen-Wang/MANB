clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'functions'));
%% ====== Parameter Settings ======
rho = 1000;       % Medium density (e.g., water: 1000 kg/m^3)
c   = 1000;       % Sound speed (m/s)
d   = 100;          % Reference period length (for normalization)
k0  = 2*pi / d;   % Reference wavenumber

% Different sequence lengths
TLn_fib = 8;
TLn_tm  = 5;
TLn_dp  =5;

% Cross-sectional area for each type of layer
S_A = 1;
S_B = 2;

%% ====== Generate Numeric Sequences ======
fib_seq = fibonacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end




























onacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end



























onacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end




























onacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end




























onacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end



























onacci_numeric(TLn_fib, S_A, S_B);
tm_seq  = thuemorse_numeric(TLn_tm, S_A, S_B);
dp_seq  = doubleperiod_numeric(TLn_dp, S_A, S_B);

% Calculate dz for each sequence (to make total length equal)
dz1 = d / length(fib_seq);
dz2 = d / length(tm_seq);
dz3 = d / length(dp_seq);

%% ====== Set k/k0 Range ======
k_over_k0_max = 16.1;  % Maximum value of k/k0
k_range = linspace(0, k_over_k0_max * k0, 2000);  % k range

% Initialize storage arrays
TL_fib = zeros(size(k_range));
TL_tm  = zeros(size(k_range));
TL_dp  = zeros(size(k_range));

%% ====== Calculate TL ======
fprintf('Calculating TL for Fibonacci sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_fib(i), ~] = calcTL_and_t_from_seq_array(fib_seq, k, dz1, rho, c);
end

fprintf('Calculating TL for Thue-Morse sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_tm(i), ~] = calcTL_and_t_from_seq_array(tm_seq, k, dz2, rho, c);
end

fprintf('Calculating TL for Double Period sequence...\n');
for i = 1:length(k_range)
    k = k_range(i);
    [TL_dp(i), ~] = calcTL_and_t_from_seq_array(dp_seq, k, dz3, rho, c);
end

k_over_k0_max1= 16.1;  % Maximum value of k/k0
k_range1 = linspace(0, k_over_k0_max, 2000);  % k range

AAA=[k_range1',TL_fib',TL_tm',TL_dp'];

%% ====== Plotting: TL Comparison of Three Sequences ======
figure('Position', [300, 200, 1000, 600]);

% Calculate k/k0 for horizontal axis
k_over_k0 = k_range / k0;

% Plot three curves
plot(k_over_k0, TL_fib, '-', 'LineWidth', 2, 'DisplayName', 'Fibonacci'); hold on;
plot(k_over_k0, TL_tm, '-', 'LineWidth', 2, 'DisplayName', 'Thue-Morse');
plot(k_over_k0, TL_dp, '-', 'LineWidth', 2, 'DisplayName', 'Double Period');

% Labels and title
xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('TL (dB)', 'FontSize', 24);
title('Transmission Loss Comparison of Different Sequences', 'FontSize', 26);

% Grid and legend
grid on;
legend('Location', 'best', 'FontSize', 18);

% Set axis font size
set(gca, 'FontSize', 18);

% Set y-axis range
TL_max = max([max(TL_fib), max(TL_tm), max(TL_dp)]);
ylim([0, ceil(TL_max) + 1]);

%% ====== Save as EPS ======
% Save to output directory
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
moire_dir = fileparts(figure_code_dir);  % moire folder (if exists)
output_dir = fullfile(moire_dir, 'output');
if ~exist(output_dir, 'dir')
    output_dir = figure_code_dir;  % If directory doesn't exist, save to figure_code folder
end
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences.png'));

fprintf('Calculation completed!\n');
fprintf('Fibonacci sequence length: %d\n', length(fib_seq));
fprintf('Thue-Morse sequence length: %d\n', length(tm_seq));
fprintf('Double Period sequence length: %d\n', length(dp_seq));

%% ====== Narrow Band Detection and Sensitivity Calculation ======
fprintf('\n=== Starting Narrow Band Detection and Sensitivity Calculation ===\n\n');

% Narrow band detection parameters
k_min0 = 0;
k_max0 = k_over_k0_max * k0;
num_k_coarse = 2001;
num_k_fine = 100;

% Detection thresholds
slope_threshold = pi;        % Absolute value threshold for slope
width_threshold = 0.5;        % Δk/k0 threshold
amplitude_threshold = 1;     % Narrow band amplitude threshold (max - min)

% Multiple refinement scanning parameters
max_refinement_iter =15;     % Maximum number of refinement scanning iterations
refinement_margin = 0.3;     % Boundary expansion ratio per iteration
refinement_points_factor = 1.5; % Points multiplication factor per iteration

% Sensitivity calculation parameters
n_T =1;               % Number of periods (for transmission calculation)
delta_k = 1e-4;     % Frequency variation range (for sensitivity calculation)
num_k_sens = 501;        % Number of frequency scanning points
delta=1e-3;
% Store results
sequences = {fib_seq, tm_seq, dp_seq};
seq_names = {'Fibonacci', 'Thue-Morse', 'Double Period'};
seq_dzs = [dz1, dz2, dz3];
final_results = cell(3, 1);
refinement_histories = cell(3, 1);
candidates_all = cell(3, 1);
k_coarse_all = cell(3, 1);
real_coarse_all = cell(3, 1);
k_min = zeros(3, 1);
k_max = zeros(3, 1);
sensitivity_results = cell(3, 1);
sensitivity_delta_k_all = cell(3, 1);
sensitivity_intensity_all = cell(3, 1);
k_final_all = cell(3, 1);
real_final_all = cell(3, 1);

%% ====== Find Narrow Bands for Three Sequences ======
for seq_idx = 1:3
    fprintf('\n========== %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;  % Use same reference length
    
    try
        [final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
            find_narrow_band_from_seq(...
            S_seq, d_seq, k0, ...
            k_min0, k_max0, ...
            num_k_coarse, num_k_fine, ...
            slope_threshold, width_threshold, amplitude_threshold, ...
            max_refinement_iter, refinement_margin, refinement_points_factor, ...
            rho, c, ...
            true);  % verbose = true
        
        final_results{seq_idx} = final_result;
        refinement_histories{seq_idx} = refinement_history;
        candidates_all{seq_idx} = candidates;
        k_coarse_all{seq_idx} = k_coarse;
        real_coarse_all{seq_idx} = real_coarse;
        
        % Calculate center frequency
        k_min(seq_idx) = (final_result(1)+final_result(2))/2 ;
        k_max(seq_idx)=k_min(seq_idx)+delta;
        
        % Perform final refinement scan for plotting
        kf_final_min = max(final_result(1) - 0.1*(final_result(2)-final_result(1)), k_min0);
        kf_final_max = min(final_result(2) + 0.1*(final_result(2)-final_result(1)), k_max0);
        [real_final, ~, k_final] = run_scan_from_seq(kf_final_min, kf_final_max, num_k_fine*2, S_seq, d_seq, rho, c);
        k_final_all{seq_idx} = k_final;
        real_final_all{seq_idx} = real_final;
        
         fprintf('\n%s Narrow Band Results:\n', seq_names{seq_idx});
        fprintf('  k_min = %.12e (%.6f k0)\n', final_result(1), final_result(1)/k0);
        fprintf('  k_max = %.12e (%.6f k0)\n', final_result(2), final_result(2)/k0);
        
        fprintf('  Δk/k0 = %.12e\n', final_result(3));
        
    catch ME
        fprintf('Warning: %s sequence narrow band detection failed: %s\n', seq_names{seq_idx}, ME.message);
        final_results{seq_idx} = [];
    end
end

%% ====== Calculate Sensitivity ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    fprintf('\n========== Calculating Sensitivity for %s Sequence ==========\n', seq_names{seq_idx});
    S_seq = sequences{seq_idx};
    d_seq = d;
    
    try
        [sensitivity_result, delta_k_values, intensity2] = ...
            calculate_sensitivity_from_seq(...
            k_min(seq_idx), S_seq, d_seq, rho, c, n_T, k_max(seq_idx), num_k_sens, ...
            [], false);  % Don't save separate sensitivity figure
        sensitivity_results{seq_idx} = sensitivity_result;
        sensitivity_delta_k_all{seq_idx} = delta_k_values;
        sensitivity_intensity_all{seq_idx} = intensity2;
        
        fprintf('%s Sensitivity Results:\n', seq_names{seq_idx});
        fprintf('  Maximum sensitivity: %.6e\n', sensitivity_result.sensitivity_max);
        fprintf('  Peak position: δk = %.6e\n', sensitivity_result.sensitivity_peak_delta_k);
        
    catch ME
        fprintf('Warning: %s sequence sensitivity calculation failed: %s\n', seq_names{seq_idx}, ME.message);
        sensitivity_results{seq_idx} = [];
    end
end

%% ====== Plotting: Band Structure and Sensitivity Comparison (Combined into One Large Figure) ======
figure('Position', [100, 100, 1600, 1000]);

% Subplot label letters
band_labels = {'(a)', '(b)', '(c)'};
sens_labels = {'(d)', '(e)', '(f)'};

% ====== First Row: Band Structure and Narrow Band Marking ======
for seq_idx = 1:3
    if isempty(final_results{seq_idx})
        continue;
    end
    
    % Calculate main plot position (2 rows, 3 columns, first row)
    subplot(2, 3, seq_idx);
    
    % Plot band structure, add legend label
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};
    plot(k_coarse/k0, real_coarse(1,:), 'b-', 'LineWidth', 2.25, 'DisplayName', sprintf('%s', seq_names{seq_idx})); hold on;
    plot(k_coarse/k0, real_coarse(2,:), 'r-', 'LineWidth', 2.25, 'HandleVisibility', 'off');
    
    % Mark candidate narrow bands
    candidates = candidates_all{seq_idx};
    y_data_min = min(real_coarse(:));
    y_data_max = max(real_coarse(:));
    
    % Sort candidates by average slope
    if ~isempty(candidates)
        [~, sort_idx] = sort(candidates(:,4), 'descend');
        best_candidate_idx = sort_idx(1);
        
        for i = 1:size(candidates,1)
            color = [0.9, 0.9, 0.9];
            if i == best_candidate_idx
                color = [0.8, 1, 0.8];
            end
            fill([candidates(i,1)/k0, candidates(i,2)/k0, candidates(i,2)/k0, candidates(i,1)/k0], ...
                 [y_data_min, y_data_min, y_data_max, y_data_max], ...
                 color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    % Mark iteration history
    refinement_history = refinement_histories{seq_idx};
    if ~isempty(refinement_history) && size(refinement_history, 1) > 1
        colors = lines(size(refinement_history,1));
        for i = 1:size(refinement_history,1)
            plot([refinement_history(i,1)/k0, refinement_history(i,1)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot([refinement_history(i,2)/k0, refinement_history(i,2)/k0], [y_data_min, y_data_max], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Frame final narrow band with green box
    final_result = final_results{seq_idx};
    rect_left = final_result(1)/k0 - 0.25;
    rect_width = (final_result(2)-final_result(1))/k0 + 0.5;
    rectangle('Position', [rect_left, y_data_min, rect_width, y_data_max-y_data_min], ...
              'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');
    
    % Set axis ranges
    xlim([min(k_coarse/k0), max(k_coarse/k0)]);
    ylim([-pi, pi]);
    
    % Set y-axis ticks as multiples of pi
    yticks([-pi, 0, pi]);
    yticklabels({'$-\pi$', '$0$', '$\pi$'});
    
    % Set x-axis ticks
    xticks([0, 5, 10, 15]);
    xticklabels({'$0$', '$5$', '$10$', '$15$'});
    
    xlabel('$k/k_0$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$qd$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Add subplot label (top left)
    text(0.02, 0.98, band_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add refinement scan subplot (at top left of each main plot)
    if ~isempty(k_final_all{seq_idx})
        % Get current subplot position information, calculate absolute position of subplot
        pos = get(gca, 'Position');
        axes('Position', [pos(1)+0.01, pos(2)+pos(4)*0.55, pos(3)*0.55, pos(4)*0.40]);
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};
        plot(k_final/k0, real_final(1,:), 'b-', 'LineWidth', 2.25); hold on;
        plot(k_final/k0, real_final(2,:), 'r-', 'LineWidth', 2.25);
        xlim([min(k_final/k0), max(k_final/k0)]);
        ylim([min(real_final(:)), max(real_final(:))]);
        set(gca, 'XTickLabel', []);  % Remove x-axis labels
        set(gca, 'YTickLabel', []);  % Remove y-axis labels
        grid on;
        set(gca, 'FontSize', 11);
    end
end

% ====== Second Row: Sensitivity Comparison ======
for seq_idx = 1:3
    if isempty(sensitivity_results{seq_idx})
        continue;
    end
    
    % Calculate subplot position (2 rows, 3 columns, second row)
    subplot(2, 3, seq_idx + 3);
    
    delta_k_values = sensitivity_delta_k_all{seq_idx};
    intensity2 = sensitivity_intensity_all{seq_idx};
    sensitivity_result = sensitivity_results{seq_idx};
    
    % Refer to Moire method: first calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only plot delta_k > 0 part
    positive_idx = delta_k_values > 0;
    delta_k_positive = delta_k_values(positive_idx);
    log_intensity_positive = log_intensity(positive_idx);
    sensitivity_positive = sensitivity_result.sensitivity(positive_idx);
    
    % Left y-axis: logarithmic transmission intensity, add legend label
    yyaxis left
    plot(delta_k_positive, log_intensity_positive, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s', seq_names{seq_idx}));
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_positive, sensitivity_positive, 'r-', 'LineWidth', 2.5, ...
         'HandleVisibility', 'off');
    ylabel('Sensitivity', 'FontSize', 20);
    set(gca, 'YColor', 'r');
    
    % x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Add subplot label (top left)
    text(0.02, 0.98, sens_labels{seq_idx}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Add legend
    if seq_idx == 1
        legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    end
    
    grid on;
    set(gca, 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
end

% Save combined figure
print('-depsc', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));
fprintf('\nCombined band structure and sensitivity comparison figure saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.eps'));
fprintf('PNG format image saved to: %s\n', fullfile(output_dir, 'Fig_TL_Sequences_Band_Sensitivity_Combined.png'));

%% ====== Unified Data Export to Current Folder/figure ======

script_dir = fileparts(mfilename('fullpath'));
figure_dir = fullfile(script_dir, 'figure');
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fprintf('\n=== Starting to Export All CSV Data ===\n');

% ----------------------------
% 1. Export TL Data
% ----------------------------
csv_file_tl = fullfile(figure_dir, 'TL_Sequences_Data.csv');

% Ensure column vectors
k_over_k0_col = k_over_k0(:);
TL_fib_col = TL_fib(:);
TL_tm_col  = TL_tm(:);
TL_dp_col  = TL_dp(:);

% Align (take maximum length, pad with NaN)
L = max([length(k_over_k0_col), length(TL_fib_col), length(TL_tm_col), length(TL_dp_col)]);
ktmp = nan(L,1); ktmp(1:length(k_over_k0_col)) = k_over_k0_col;
a = nan(L,1); a(1:length(TL_fib_col)) = TL_fib_col;
b = nan(L,1); b(1:length(TL_tm_col)) = TL_tm_col;
c = nan(L,1); c(1:length(TL_dp_col)) = TL_dp_col;

T_tl = table(ktmp, a, b, c, 'VariableNames', {'k_over_k0','TL_Fibonacci','TL_ThueMorse','TL_DoublePeriod'});
writetable(T_tl, csv_file_tl);
fprintf('Exported TL data: %s\n', csv_file_tl);


% ----------------------------
% 2. Export Band Data (Coarse Scan + Fine Scan)
% ----------------------------
seq_names = {'Fibonacci','ThueMorse','DoublePeriod'};

for seq_idx = 1:3

    % Skip failed sequences
    if isempty(k_coarse_all{seq_idx})
        fprintf('Skipping %s: No coarse scan data\n', seq_names{seq_idx});
        continue;
    end

    % ---- (A) Coarse Scan Bands ----
    k_coarse = k_coarse_all{seq_idx};
    real_coarse = real_coarse_all{seq_idx};

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
                error('Cannot identify real_coarse size (neither 2xN nor Nx2)');
            end
        end
    catch ME
        warning('Sequence %s: Failed to parse real_coarse: %s. Skipping coarse scan band export.', seq_names{seq_idx}, ME.message);
        continue;
    end

    % Align (pad to maximum length with NaN)
    kc = k_coarse(:);
    Lc = max([length(kc), length(band1), length(band2)]);
    kc_pad = nan(Lc,1); kc_pad(1:length(kc)) = kc;
    b1_pad = nan(Lc,1); b1_pad(1:length(band1)) = band1;
    b2_pad = nan(Lc,1); b2_pad(1:length(band2)) = band2;

    coarse_file = fullfile(figure_dir, sprintf('Band_%s_Coarse.csv', seq_names{seq_idx}));
    T_band_coarse = table(kc_pad./k0, b1_pad, b2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
    writetable(T_band_coarse, coarse_file);
    fprintf('Exported coarse scan bands: %s\n', coarse_file);


    % ---- (B) Fine Scan Bands (if exists) ----
    if ~isempty(k_final_all{seq_idx}) && ~isempty(real_final_all{seq_idx})
        k_final = k_final_all{seq_idx};
        real_final = real_final_all{seq_idx};

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
                    error('Cannot identify real_final size (neither 2xN nor Nx2)');
                end
            end
        catch ME
            warning('Sequence %s: Failed to parse real_final: %s. Skipping fine scan band export.', seq_names{seq_idx}, ME.message);
            continue;
        end

        % Align (pad to maximum length with NaN)
        kf = k_final(:);
        Lf = max([length(kf), length(fband1), length(fband2)]);
        kf_pad = nan(Lf,1); kf_pad(1:length(kf)) = kf;
        fb1_pad = nan(Lf,1); fb1_pad(1:length(fband1)) = fband1;
        fb2_pad = nan(Lf,1); fb2_pad(1:length(fband2)) = fband2;

        final_file = fullfile(figure_dir, sprintf('Band_%s_Final.csv', seq_names{seq_idx}));
        T_band_final = table(kf_pad./k0, fb1_pad, fb2_pad, 'VariableNames', {'k_over_k0','band1','band2'});
        writetable(T_band_final, final_file);
        fprintf('Exported fine scan bands: %s\n', final_file);
    else
        fprintf('Sequence %s: No fine scan data, skipping Final export\n', seq_names{seq_idx});
    end
end


% ----------------------------
% 3. Export Sensitivity Data (1 file per sequence)
% ----------------------------
for seq_idx = 1:3

    if isempty(sensitivity_results{seq_idx}) || isempty(sensitivity_delta_k_all{seq_idx}) || isempty(sensitivity_intensity_all{seq_idx})
        fprintf('Sequence %s: No sensitivity data or calculation failed, skipping\n', seq_names{seq_idx});
        continue;
    end

    delta_k_values = sensitivity_delta_k_all{seq_idx}(:);
    intensity2 = sensitivity_intensity_all{seq_idx}(:);
    
    % Ensure data are column vectors
    if isrow(delta_k_values)
        delta_k_values = delta_k_values';
    end
    if isrow(intensity2)
        intensity2 = intensity2';
    end
    
    % Get sensitivity data
    sens = sensitivity_results{seq_idx}.sensitivity;
    if isrow(sens)
        sens = sens';
    end
    
    % Check data length
    % After using gradient, sensitivity length is same as delta_k_values
    % If lengths don't match, perform alignment
    if length(sens) < length(delta_k_values)
        % If sensitivity length is shorter, add NaN padding
        sens = [sens; nan(length(delta_k_values) - length(sens), 1)];
    elseif length(sens) > length(delta_k_values)
        % If sensitivity length is longer, truncate
        sens = sens(1:length(delta_k_values));
    end
    
    % Refer to Moire method: calculate logarithmic normalized pressure magnitude squared
    log_intensity = log(intensity2 / min(intensity2));
    
    % Only keep delta_k > 0 part (consistent with plotting range)
    positive_idx = delta_k_values > 0;
    delta_k_values = delta_k_values(positive_idx);
    intensity2 = intensity2(positive_idx);
    sens = sens(positive_idx);
    log_intensity = log_intensity(positive_idx);
    
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

    sens_file = fullfile(figure_dir, sprintf('Sensitivity_%s.csv', seq_names{seq_idx}));
    T_sens = table(dk_pad, int_pad, log_int_pad, sens_pad, ...
        'VariableNames', {'delta_k','intensity','log_intensity','sensitivity'});
    writetable(T_sens, sens_file);
    fprintf('Exported sensitivity data: %s\n', sens_file);
end