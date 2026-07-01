function ipr = calculate_ipr_at_k0(S_seq, k0, d, rho, c)
% Calculate size-independent IPR (Inverse Participation Ratio) at k0
%
% Parameters:
%   S_seq: Cross-section sequence array (column or row vector)
%   k0: Target wavenumber
%   d: Total length
%   rho: Density
%   c: Sound speed
%
% Returns:
%   ipr: Size-independent IPR value

% Add functions directory to path
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

num_n = length(S_seq);
dz = d / num_n;

% Build transfer matrix
TM_T = eye(2);
for jj = 1:num_n
    S_curr = S_seq(jj);
    T1 = acoustic_TM(S_curr, k0, dz, rho, c);
    TM_T = T1 * TM_T;
end

% Calculate eigenvalues and eigenvectors
[V, D] = eig(TM_T);
eigval = diag(D);

% Calculate Bloch phase
revals = real(-1i * log(eigval));
imevals = imag(-1i * log(eigval));
ime_copy = imevals;
revals = mod(revals + pi, 2*pi) - pi;
re_copy = revals;

% Sort eigenvalues
if revals(2) > revals(1)
    revals(1) = re_copy(2);
    revals(2) = re_copy(1);
    % Simultaneously swap eigenvectors
    V_temp = V(:, 1);
    V(:, 1) = V(:, 2);
    V(:, 2) = V_temp;
end
revals(1) = -revals(2);

if imevals(2) > imevals(1)
    imevals(1) = ime_copy(2);
    imevals(2) = ime_copy(1);
end

% Select first eigenvector (corresponds to band 1)
psi = V(:, 1);

% Normalize
psi = psi / norm(psi);

% Calculate size-independent IPR
ipr = calculate_size_independent_ipr(psi);

end

