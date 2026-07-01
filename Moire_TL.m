clc; clear; close all;

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);  % figure_code folder
addpath(fullfile(figure_code_dir, 'functions'));  % Use local functions folder

%% ====== Parameter Settings ======
rho = 1000;       % Medium density (e.g., water: 1000 kg/m^3)
c   = 1000;       % Sound speed (m/s)
N1 = 12;          % First subsection sequence length
N2 = N1 + 1;      % Second subsection length (Moiré effect)
d  = lcm(N1, N2); % Basic period length
num_n = d * 8;    % Subdivide one period into num_n segments
dz = d / num_n;   % Thickness per segment
k0 = 2*pi / d;    % Reference wavenumber
Sa = 1; Sb = 4;   % Subsection cross-sectional area
Sa1 = 1; Sb1 = 4;
kd_range = linspace(0, 11*k0, 1001); % kd range

%% ====== Generate Moiré Sequence ======
n1_t = gener_n(num_n, N2, Sa, Sb);
n2_t = gener_n(num_n, N1, Sa1, Sb1);
n_moire = n1_t + n2_t; % Moiré superposition sequence

%% ====== Calculate TL and |t| ======
TL_moire = zeros(size(kd_range));
t_abs = zeros(size(kd_range));

for i = 1:length(kd_range)
    k = kd_range(i);
    [TL_moire(i), t_abs(i)] = calcTL_and_t_from_seq_array(n_moire, k, dz, rho, c);
end

%% ====== Narrow Band Position Parameters ======
% Obtained from narrow band detection results (refer to test3.m)
k_min_nb = 2.422789740192e-01-0.00005;  % Narrow band minimum wavenumber
k_max_nb = 2.422805299970e-01+0.00005;  % Narrow band maximum wavenumber
k_center_nb = 6.01536*k0;  % Narrow band center wavenumber

% Add 2000 points near narrow band center (+-0.5k0) for more precise plotting (add to main plot)
k_nb_fine_range = linspace(k_center_nb - 0.5*k0, k_center_nb + 0.5*k0,100);
TL_nb_fine = zeros(size(k_nb_fine_range));
t_abs_nb_fine = zeros(size(k_nb_fine_range));

fprintf('Calculating fine scan near narrow band center (2000 points, range: %.6f to %.6f)...\n', ...
        k_center_nb - 0.5*k0, k_center_nb + 0.5*k0);
for i = 1:length(k_nb_fine_range)
    k = k_nb_fine_range(i);
    [TL_nb_fine(i), t_abs_nb_fine(i)] = calcTL_and_t_from_seq_array(n_moire, k, dz, rho, c);
end
fprintf('Fine scan completed!\n');

% Calculate fine scan near narrow band (for right subplot)
k_nb_range = linspace(k_min_nb - 0.1*(k_max_nb-k_min_nb), ...
                       k_max_nb + 0.1*(k_max_nb-k_min_nb), 1000);
TL_nb = zeros(size(k_nb_range));
t_abs_nb = zeros(size(k_nb_range));

for i = 1:length(k_nb_range)
    k = k_nb_range(i);
    [TL_nb(i), t_abs_nb(i)] = calcTL_and_t_from_seq_array(n_moire, k, dz, rho, c);
end

%% ====== Calculate Data Points at Narrow Band Position ======
% Find data points corresponding to narrow band position in main plot
[~, idx_min] = min(abs(kd_range - k_min_nb));
[~, idx_max] = min(abs(kd_range - k_max_nb));
[~, idx_center] = min(abs(kd_range - k_center_nb));

% Calculate TL and |t| values at narrow band position
TL_nb_min = TL_moire(idx_min);
TL_nb_max = TL_moire(idx_max);
TL_nb_center = TL_moire(idx_center);
t_nb_min = t_abs(idx_min);
t_nb_max = t_abs(idx_max);
t_nb_center = t_abs(idx_center);

%% ====== Plotting: Main Plot + Narrow Band Zoom ======
figure('Position',[300 200 1600 600]);

% ====== Left Plot: Main Plot (Full Range) ======
subplot(1, 2, 1);
% Adjust subplot position to reduce blank spacing
pos1 = get(gca, 'Position');
set(gca, 'Position', [0.08, 0.15, 0.42, 0.75]);

% Left y-axis (TL)
yyaxis left
scatter(kd_range/k0, TL_moire, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6); hold on;
% Add fine scan data near narrow band center (2000 points, add to main plot)
scatter(k_nb_fine_range/k0, TL_nb_fine, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.9);
% Add subplot scan data to main plot
scatter(k_nb_range/k0, TL_nb, 25, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
ylabel('TL (dB)', 'FontSize', 24);

% Left y-axis only shows 3 values (min, center, max), and as integers
TL_max = ceil(max([max(TL_moire), max(TL_nb_fine), max(TL_nb)]));           % Round up maximum value
TL_min = 0;
TL_center = round(TL_max / 2);  % Center value also rounded
yticks([TL_min, TL_center, TL_max]);
yticklabels({num2str(TL_min), num2str(TL_center), num2str(TL_max)});  % Ensure displayed as integers
ylim([0 TL_max]);

% Right y-axis (|t|)
yyaxis right
scatter(kd_range/k0, t_abs, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6); hold on;
% Add fine scan data near narrow band center (2000 points, add to main plot)
scatter(k_nb_fine_range/k0, t_abs_nb_fine, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.9);
% Add subplot scan data to main plot
scatter(k_nb_range/k0, t_abs_nb, 25, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
ylabel('|t| ', 'FontSize', 28);
ylim([0, 1]);
yticks([0, 0.5, 1]);  % Only show 3 values: left, center, right

% Frame narrow band region with green box
y_min_rect = 0;
y_max_rect = max([TL_max, 1]);  % Take maximum of TL and |t|
rect_left = k_min_nb/k0 - 0.25;  % Left boundary extends 0.25*k0
rect_width = (k_max_nb - k_min_nb)/k0 + 0.5;  % Width increases by 0.5*k0
rectangle('Position', [rect_left, y_min_rect, rect_width, y_max_rect - y_min_rect], ...
          'EdgeColor', 'g', 'LineWidth', 2.5, 'LineStyle', '-', 'HandleVisibility', 'off');

% x-axis - only show 3 values (min, center, max)
x_min = 0;
x_max = max(kd_range/k0);
x_center = x_max / 2;
xticks([x_min, x_center, x_max]);
xlabel('$k/k_0$', 'Interpreter','latex', 'FontSize', 24);

% Add subplot label (a) at top left
text(0.02, 0.98, '(a)', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
     'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Grid and legend
grid on;
legend({'TL (dB)', '|t|'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');

% Set font size for entire plot
set(gca, 'FontSize', 18);
set(gca, 'TickLabelInterpreter', 'latex');

% ====== Right Plot: Narrow Band Region Zoom (Dual y-axis) ======
subplot(1, 2, 2);
% Adjust subplot position to reduce blank spacing
pos2 = get(gca, 'Position');
set(gca, 'Position', [0.55, 0.15, 0.42, 0.75]);

% Left y-axis (TL)
yyaxis left
scatter(k_nb_range/k0, TL_nb, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.8); hold on;
ylabel('TL (dB)', 'FontSize', 24);
TL_nb_max_val = ceil(max(TL_nb)*1.1);  % Round up
TL_nb_min_val = 0;
TL_nb_center_val = round(TL_nb_max_val / 2);  % Center value also rounded
ylim([TL_nb_min_val, TL_nb_max_val]);
yticks([TL_nb_min_val, TL_nb_center_val, TL_nb_max_val]);  % Only show 3 values
yticklabels({num2str(TL_nb_min_val), num2str(TL_nb_center_val), num2str(TL_nb_max_val)});  % Ensure displayed as integers

% Right y-axis (|t|)
yyaxis right
scatter(k_nb_range/k0, t_abs_nb, 30, 'r', 'filled', 'MarkerFaceAlpha', 0.8); hold on;
ylabel('|t| ', 'FontSize', 28);
ylim([0, 1]);
yticks([0, 0.5, 1]);  % Only show 3 values: left, center, right

% Mark narrow band boundaries (green dashed lines)
yyaxis left
plot([k_min_nb/k0, k_min_nb/k0], [0, max(TL_nb)*1.1], 'g--', 'LineWidth', 2, ...
     'HandleVisibility', 'off');
plot([k_max_nb/k0, k_max_nb/k0], [0, max(TL_nb)*1.1], 'g--', 'LineWidth', 2, ...
     'HandleVisibility', 'off');

yyaxis right
plot([k_min_nb/k0, k_min_nb/k0], [0, 1], 'g--', 'LineWidth', 2, ...
     'HandleVisibility', 'off');
plot([k_max_nb/k0, k_max_nb/k0], [0, 1], 'g--', 'LineWidth', 2, ...
     'HandleVisibility', 'off');

% x-axis - show narrow band boundary coordinate values and center value (3 values)
x_nb_min = k_min_nb/k0;
x_nb_max = k_max_nb/k0;
x_nb_center = k_center_nb/k0;
xticks([x_nb_min, x_nb_center, x_nb_max]);  % Show coordinate values of green lines
xlabel('$k/k_0$', 'Interpreter','latex', 'FontSize', 24);

% Add subplot label (b) at top left
text(0.02, 0.98, '(b)', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
     'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Grid and legend
grid on;
legend({'TL (dB)', '|t|'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');

% Set font size for entire plot
set(gca, 'FontSize', 18);
set(gca, 'TickLabelInterpreter', 'latex');

%% ====== Save Plot Data to figure Folder (CSV Format) ======
% Create or check figure folder
output_dir = fullfile(script_dir, 'figure');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created figure folder: %s\n', output_dir);
end

% 1. Save main plot data (full range)
k_over_k0_main = kd_range(:) / k0;
TL_main = TL_moire(:);
t_abs_main = t_abs(:);
T_main = table(k_over_k0_main, TL_main, t_abs_main, ...
    'VariableNames', {'k_over_k0', 'TL_dB', 't_abs'});
csv_file_main = fullfile(output_dir, 'Moire_TL_main.csv');
writetable(T_main, csv_file_main);
fprintf('Main plot data saved to: %s\n', csv_file_main);

% 2. Save narrow band fine scan data (for main plot overlay)
k_over_k0_nb_fine = k_nb_fine_range(:) / k0;
TL_nb_fine_col = TL_nb_fine(:);
t_abs_nb_fine_col = t_abs_nb_fine(:);
T_nb_fine = table(k_over_k0_nb_fine, TL_nb_fine_col, t_abs_nb_fine_col, ...
    'VariableNames', {'k_over_k0', 'TL_dB', 't_abs'});
csv_file_nb_fine = fullfile(output_dir, 'Moire_TL_nb_fine.csv');
writetable(T_nb_fine, csv_file_nb_fine);
fprintf('Narrow band fine scan data saved to: %s\n', csv_file_nb_fine);

% 3. Save narrow band region scan data (for right subplot)
k_over_k0_nb = k_nb_range(:) / k0;
TL_nb_col = TL_nb(:);
t_abs_nb_col = t_abs_nb(:);
T_nb = table(k_over_k0_nb, TL_nb_col, t_abs_nb_col, ...
    'VariableNames', {'k_over_k0', 'TL_dB', 't_abs'});
csv_file_nb = fullfile(output_dir, 'Moire_TL_nb.csv');
writetable(T_nb, csv_file_nb);
fprintf('Narrow band region scan data saved to: %s\n', csv_file_nb);
