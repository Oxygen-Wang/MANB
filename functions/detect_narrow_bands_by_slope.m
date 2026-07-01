function candidates = detect_narrow_bands_by_slope(k_pro, real_eigval_pro, k0, slope_threshold, width_threshold, amplitude_threshold)
% Detect narrow bands using slope method on given data
%
% Parameters:
%   k_pro: Wavenumber array
%   real_eigval_pro: Real part of eigenvalues (2 x num_k)
%   k0: Normalization wavenumber
%   slope_threshold: Slope threshold
%   width_threshold: Width threshold
%   amplitude_threshold: Amplitude threshold
%
% Returns:
%   candidates: Candidate narrow band array [k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude, index]

q = real_eigval_pro(1,:);
dq = diff(q);
dk = diff(k_pro);
slope = abs(dq ./ dk);

% Find high slope regions
high_slope_mask = slope > slope_threshold;
high_slope_idx = find(high_slope_mask);

if isempty(high_slope_idx)
    candidates = [];
    return;
end

% Group contiguous high slope regions
slope_groups = contiguous_groups(high_slope_idx);

candidates = [];
for i = 1:size(slope_groups, 1)
    idx_start = max(slope_groups(i, 1), 1);
    idx_end = min(slope_groups(i, 2) + 1, length(k_pro));
    
    if idx_end > idx_start
        k_min = k_pro(idx_start);
        k_max = k_pro(idx_end);
        delta_k_norm = (k_max - k_min) / k0;
        
        % Calculate average slope in this interval
        local_slopes = slope(idx_start:idx_end-1);
        avg_slope = mean(local_slopes);
        max_slope = max(local_slopes);
        
        % Calculate amplitude of band values within narrow band (max - min)
        band_values = q(idx_start:idx_end);
        amplitude = max(band_values) - min(band_values);
        
        % Only keep candidates satisfying width and amplitude conditions
        if delta_k_norm < width_threshold && amplitude >= amplitude_threshold
            candidates = [candidates; k_min, k_max, delta_k_norm, avg_slope, max_slope, amplitude, i];
        end
    end
end

end













