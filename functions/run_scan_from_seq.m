function [real_eigval_pro, imag_eigval_pro, k_pro] = run_scan_from_seq(k_min, k_max, num_k, S_seq, d, rho, c)
% Perform band structure scan from sequence array
%
% Parameters:
%   k_min, k_max: Scanning range
%   num_k: Number of scanning points
%   S_seq: Cross-section sequence array (column or row vector)
%   d: Total length (for normalization)
%   rho: Density
%   c: Sound speed
%
% Returns:
%   real_eigval_pro: Real part of eigenvalues (2 x num_k)
%   imag_eigval_pro: Imaginary part of eigenvalues (2 x num_k)
%   k_pro: Wavenumber array

k_pro = linspace(k_min, k_max, num_k);
real_eigval_pro = zeros(2, num_k);
imag_eigval_pro = zeros(2, num_k);

num_n = length(S_seq);
dz = d / num_n;

for ii = 1:num_k
    k = k_pro(ii);
    TM_T = eye(2);
    for jj = 1:num_n
        S_curr = S_seq(jj);
        T1 = acoustic_TM(S_curr, k, dz, rho, c);
        TM_T = T1 * TM_T;
    end
    eigval = eig(TM_T);
    revals = real(-1i * log(eigval));
    imevals = imag(-1i * log(eigval));
    ime_copy = imevals;
    revals = mod(revals + pi, 2*pi) - pi;
    re_copy = revals;
    if revals(2) > revals(1)
        revals(1) = re_copy(2);
        revals(2) = re_copy(1);
    end
    revals(1) = -revals(2);
    if imevals(2) > imevals(1)
        imevals(1) = ime_copy(2);
        imevals(2) = ime_copy(1);
    end
    real_eigval_pro(:,ii) = revals;
    imag_eigval_pro(:,ii) = imevals;
end

end

