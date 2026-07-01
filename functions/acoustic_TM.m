function Mj = acoustic_TM(S_curr, k, dz, rho, c)
% Calculate acoustic transfer matrix
%
% Parameters:
%   S_curr: Current cross-sectional area
%   k: Wavenumber
%   dz: Length increment
%   rho: Density
%   c: Sound speed
%
% Returns:
%   Mj: 2x2 transfer matrix
%
% Example:
%   Mj = acoustic_TM(0.01, 1.0, 0.1, 1000, 1500);

Zc = rho * c / S_curr;

Mj = [cos(k*dz), -1i*Zc*sin(k*dz);
      -1i*sin(k*dz)/Zc, cos(k*dz)];

end

