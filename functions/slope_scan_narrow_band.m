function [k_peak, S_max, k_vec, S_vec, I_vec] = slope_scan_narrow_band(...
    k_min, k_max, n_pro, d, num_n, rho, c, num_k)

dz = d / num_n;

k_vec = linspace(k_min, k_max, num_k);

I_vec = zeros(num_k,1);

for i = 1:num_k
    k = k_vec(i);

    state = get_ini(k, num_n, n_pro, d, rho, c);

    for j = 1:length(n_pro)
        TM = acoustic_TM(n_pro(j), k, dz, rho, c);
        state = TM * state;
    end

    I_vec(i) = abs(state(1))^2;
end

% =========================
% slope (no abs first)
% =========================
logI = log(I_vec + eps);

S_raw = gradient(logI, k_vec);

% smoothing (VERY important)
S_smooth = movmean(S_raw, max(3, round(num_k*0.03)));

% physical sensitivity definition
S_vec = abs(S_smooth);

% peak
[S_max, idx] = max(S_vec);
k_peak = k_vec(idx);

end