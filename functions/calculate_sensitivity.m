function [sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity(...
    k_min, n_pro, d, num_n, rho, c, n_T, k_max, num_k_sens, output_dir, verbose)
% Calculate transmission sensitivity near narrow band center frequency
%
% Parameters:
%   k_center: Narrow band center frequency
%   n_pro: Cross-section sequence (one superperiod)
%   d: Superperiod length
%   num_n: Number of sampling points
%   rho: Density
%   c: Sound speed
%   n_T: Number of periods (for constructing long structure)
%   delta_k: Frequency variation range
%   num_k_sens: Number of frequency scanning points (default 51)
%   output_dir: Output directory (for saving figures, optional)
%   verbose: Whether to display detailed information (default true)
%
% Returns:
%   sensitivity_result: Structure containing sensitivity-related results
%       - sensitivity: Sensitivity array
%       - sensitivity_max: Maximum sensitivity
%       - sensitivity_peak_idx: Peak position index
%       - sensitivity_peak_delta_k: delta_k corresponding to peak position
%       - intensity2: Transmission intensity after symmetry processing
%   delta_k_values: Frequency offset value array
%   intensity2: Transmission intensity after symmetry processing
%
% Example:
%   [sens_result, delta_k_vals, intensity] = calculate_sensitivity(...
%       k_center, n_pro, d, num_n, 1000, 1000, 32, 100e-5, 51, './output', true);

if nargin < 10
    output_dir = [];
end

if nargin < 11
    verbose = true;
end

if verbose
    fprintf('Starting sensitivity calculation...\n');
    fprintf('Narrow band center frequency k_center = %.12e\n', k_min);
end

% % Calculate transmission spectrum near center frequency
% k_sens = linspace(k_min, k_max, num_k_sens);
delta1 =0*1e-1;
delta2 =(k_max-k_min)*1e-1;

num_fine = round(num_k_sens * 0.1);
num_coarse = num_k_sens - num_fine;

% 粗扫描（全局）
k_coarse = linspace(k_min, k_max, num_coarse);
% 细扫描（零附近）
k_fine = linspace(k_min-delta2, k_min + delta2, num_fine);
% 合并 + 去重 + 排序
k_sens = unique([k_coarse, k_fine]);
k_sens = sort(k_sens);
num_k_sens = length(k_sens); 



% Construct long structure (n_pro as modulation)
n_long = zeros(num_n * n_T, 1);
for ii = 1:n_T
    n_long((ii-1)*num_n+1:ii*num_n) = n_pro;
end
% Get initial conditions (eigenmode)
ini = get_ini(k_min, num_n, n_pro, d, rho, c);
% Calculate transmission intensity
intensity = zeros(num_k_sens, 1);
dz = d / num_n;
for ii = 1:num_k_sens
    k = k_sens(ii);
    state = ini;
    
    % Calculate transmission through entire long structure
    for jj = 1:length(n_long)
        S_curr = n_long(jj);
        TM_seg = acoustic_TM(S_curr, k, dz, rho, c);
        state = TM_seg * state;
    end
    
    % Calculate transmission intensity (according to theory, only use pressure component state(1))
    intensity(ii) = abs(state(1))^2;
end
intensity2 = intensity;
% k-grid relative to center
delta_k_values = k_sens - k_min;
log_intensity = log(intensity2 + eps);
dlogI_dk = gradient(log_intensity, delta_k_values);
sensitivity = dlogI_dk ;  % Multiply by total length
% Output sensitivity results
sensitivity_max = max(abs(sensitivity));
sensitivity_peak_idx = find(abs(sensitivity) == sensitivity_max, 1);
sensitivity_peak_delta_k = delta_k_values(sensitivity_peak_idx);
% Build return structure
sensitivity_result = struct();
sensitivity_result.sensitivity = sensitivity;
sensitivity_result.sensitivity_max = sensitivity_max;
sensitivity_result.sensitivity_peak_idx = sensitivity_peak_idx;
sensitivity_result.sensitivity_peak_delta_k = sensitivity_peak_delta_k;
sensitivity_result.intensity2 = intensity2;
if verbose
    fprintf('Maximum sensitivity: %.6e\n', sensitivity_max);
    fprintf('Peak position: δk = %.6e\n', sensitivity_peak_delta_k);
end
% If output directory is provided, plot and save figure
if ~isempty(output_dir)
    % Combined figure with dual y-axis: Transmission intensity and Sensitivity
    figure('Position', [100, 100, 800, 600]);
    
    % Left y-axis: Logarithmic transmission intensity
    yyaxis left
  plot(delta_k_values, log(intensity2 + eps), 'b-', 'LineWidth', 2.5);
ylabel('$\log I$', 'Interpreter', 'latex', 'FontSize', 24);
set(gca, 'YColor', 'b');

yyaxis right
plot(delta_k_values, sensitivity, 'r-', 'LineWidth', 2.5);
ylabel('Sensitivity', 'FontSize', 24);
set(gca, 'YColor', 'r');

% tight limits (no margin)
ylim_left = [min(log(intensity2 + eps)), max(log(intensity2 + eps))];
ylim_right = [min(sensitivity), max(sensitivity)];

yyaxis left;  ylim(ylim_left);
yyaxis right; ylim(ylim_right);

xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 24);

xlim([min(delta_k_values), max(delta_k_values)]);

grid on;
set(gca, 'FontSize', 20);
set(gca, 'TickLabelInterpreter', 'latex');

print('-depsc', fullfile(output_dir, 'Fig_Moire_Sensitivity_Combined.eps'));
print('-dpng', '-r300', fullfile(output_dir, 'Fig_Moire_Sensitivity_Combined.png'));
end
end

