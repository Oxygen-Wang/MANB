function [TL, t_abs] = calcTL_and_t_from_seq_array(S_array, k, dz, rho, c)
% Calculate transmission loss TL and complex amplitude transmission coefficient |t|
%
% Parameters:
%   S_array: Waveguide cross-section sequence array (column or row vector)
%   k: Wavenumber
%   dz: Length increment
%   rho: Density (kg/m^3)
%   c: Sound speed (m/s)
%
% Returns:
%   TL: Transmission loss (dB)
%   t_abs: Magnitude of complex amplitude transmission coefficient |t|
%
% Description:
%   Corresponds to paper formulas Eq.(t_general) - Eq.(TL)
%   Calculate transfer matrix for given cross-section sequence, then solve for transmission coefficient and transmission loss
%
% Example:
%   [TL, t_abs] = calcTL_and_t_from_seq_array(n_moire, k, dz, rho, c);

% Characteristic impedances at input and output ends
S_in = S_array(1);
S_out = S_array(end);
Z0 = rho * c / S_in;   % Input end characteristic impedance
ZN = rho * c / S_out;  % Output end characteristic impedance
ZN = Z0 ;  % Output end characteristic impedance
% Calculate total transfer matrix
N = length(S_array);
M = eye(2);
for i = 1:N
    M_i = acoustic_TM(S_array(i), k, dz, rho, c);
    M = M_i * M;
end

% Extract transfer matrix elements
T11 = M(1,1); T12 = M(1,2);
T21 = M(2,1); T22 = M(2,2);

% Eq.(transmission_coeff) - Calculate transmission coefficient
% According to paper formula Eq.(transmission_coeff): t = 2/(T11 + T12/ZN + Z0*T21 + T22)
% Formula: t = 2/(T_{11} + T_{12}/Z_N + Z_0 T_{21} + T_{22})
t = 2 / (T11 + T12/ZN + Z0*T21 + T22);
t_abs = abs(t);

% Eq.(TL) - Calculate transmission loss
TL = 20 * log10(1/abs(t)) + 10 * log10(ZN/Z0);



end












