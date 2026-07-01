function [sensitivity_result, delta_k_values, intensity2] = calculate_sensitivity_from_seq(...
    k_min, S_seq, d, rho, c, n_T, k_max, num_k_sens, output_dir, verbose)
% Calculate transmission sensitivity near narrow band center frequency from sequence array
%
% Parameters:
%   k_center: Narrow band center frequency
%   S_seq: Cross-section sequence array
%   d: Total length
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

if nargin < 9
    output_dir = [];
end

if nargin < 10
    verbose = true;
end

if verbose
    fprintf('Starting sensitivity calculation...\n');
    fprintf('Narrow band center frequency k_center = %.12e\n', k_min);
end

% Calculate transmission spectrum near center frequency (from k_center to k_center + delta_k)
% k_min_sens = k_center;
% k_max_sens = k_center + delta_k;
 k_sens = linspace(k_min, k_max, num_k_sens);
% delta2 =(k_max-k_min)*1e-5;
% num_fine = round(num_k_sens * 0.5);
% num_coarse = num_k_sens - num_fine;
% % 粗扫描（全局）
% k_coarse = linspace(k_min, k_max, num_coarse);
% % 细扫描（零附近）
% k_fine = linspace(k_min, k_min + delta2, num_fine);
% % 合并 + 去重 + 排序
% k_sens = unique([k_coarse, k_fine]);
% k_sens = sort(k_sens);
% num_k_sens = length(k_sens); 




% Construct long structure (S_seq as modulation)
num_n = length(S_seq);
n_long = zeros(num_n * n_T, 1);
for ii = 1:n_T
    n_long((ii-1)*num_n+1:ii*num_n) = S_seq;
end

% Get initial conditions (eigenmode)
dz = d / num_n;
TM_T = eye(2);
for jj = 1:num_n
    S_curr = S_seq(jj);
    T1 = acoustic_TM(S_curr, k_min, dz, rho, c);
    TM_T = T1 * TM_T;
end
[eigvec0,~]=eig(TM_T);
ini=eigvec0(:,1);

% Calculate transmission intensity
intensity = zeros(num_k_sens, 1);

for ii = 1:num_k_sens
    k_max = k_sens(ii);
    state = ini;
    
    % Calculate transmission through entire long structure
    for jj = 1:length(n_long)
        S_curr = n_long(jj);
        TM_seg = acoustic_TM(S_curr, k_max, dz, rho, c);
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
sensitivity = dlogI_dk ;

% Output sensitivity results (refer to test4.m)
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
    plot(delta_k_values, log_intensity_norm, 'b-', 'LineWidth', 2.5);
    ylabel('$\log |p|^2$', 'Interpreter', 'latex', 'FontSize', 24);
    set(gca, 'YColor', 'b');
    
    % Right y-axis: Sensitivity
    yyaxis right
    plot(delta_k_values, sensitivity, 'r-', 'LineWidth', 2.5);
    ylabel('Sensitivity', 'FontSize', 24);
    set(gca, 'YColor', 'r');
    
    % Common x-axis label
    xlabel('$\delta k$', 'Interpreter', 'latex', 'FontSize', 24);
    
    grid on;
    set(gca, 'FontSize', 20);
    set(gca, 'TickLabelInterpreter', 'latex');
    
    print('-depsc', fullfile(output_dir, 'Fig_Sensitivity_Combined.eps'));
    if verbose
        fprintf('Combined sensitivity figure saved to: %s\n', fullfile(output_dir, 'Fig_Sensitivity_Combined.eps'));
    end
end

end

