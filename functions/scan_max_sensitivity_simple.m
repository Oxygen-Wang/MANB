function [k_peak, S_max, k_vec, S_vec, I_vec] = scan_max_sensitivity_simple( ...
    k_min, k_max, n_pro, d, num_n, rho, c, n_T, num_k, verbose)
% 直接扫描区间内灵敏度最大点
%
% 输出:
%   k_peak : 最大灵敏度对应的 k
%   S_max  : 最大灵敏度
%   k_vec  : 扫描k
%   S_vec  : 灵敏度
%   I_vec  : 透射强度

if nargin < 10
    verbose = true;
end

dz = d / num_n;
k_vec = linspace(k_min, k_max, num_k).';
I_vec = zeros(num_k, 1);

% 构建长结构
n_long = repmat(n_pro(:), n_T, 1);

for ii = 1:num_k
    k = k_vec(ii);

    state = get_ini(k, num_n, n_pro, d, rho, c);

    for jj = 1:length(n_long)
        TM_seg = acoustic_TM(n_long(jj), k, dz, rho, c);
        state = TM_seg * state;
    end

    I_vec(ii) = abs(state(1))^2;
end

% 灵敏度定义
logI = log(I_vec + eps);
S_raw = gradient(logI, k_vec);

% 轻微平滑，避免数值尖峰
win = max(3, round(0.03 * num_k));
S_vec = abs(movmean(S_raw, win));

% 最大值
[S_max, idx] = max(S_vec);
k_peak = k_vec(idx);

if verbose
    fprintf('Sensitivity scan finished.\n');
    fprintf('Peak k = %.12e\n', k_peak);
    fprintf('Max sensitivity = %.6e\n', S_max);
end

end