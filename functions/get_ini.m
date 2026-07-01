function ini = get_ini(k, num_n, n_pro, d, rho, c)
% Get initial conditions for eigenmode
%
% Parameters:
%   k: Wavenumber
%   num_n: Number of sampling points
%   n_pro: Cross-sectional area sequence
%   d: Superperiod length
%   rho: Density
%   c: Sound speed
%
% Returns:
%   ini: Normalized eigenvector (propagating mode)
%
% Example:
%   ini = get_ini(k0, num_n, n_pro, d, 1000, 1000);

% Calculate transfer matrix for one period
dz = d / num_n;
TM_T = eye(2);
for jj = 1:num_n
    S_curr = n_pro(jj);
    T1 = acoustic_TM(S_curr, k, dz, rho, c);
    TM_T = T1 * TM_T;
end

% Calculate eigenvalues and eigenvectors
[eigvec0,eigval0]=eig(TM_T);
ini=eigvec0(:,1);
end





