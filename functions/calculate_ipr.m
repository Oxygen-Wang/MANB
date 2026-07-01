function ipr = calculate_ipr(psi)
% Calculate IPR (Inverse Participation Ratio) for a single eigenstate
%
% Parameters:
%   psi: Normalized eigenvector (column vector)
%
% Returns:
%   ipr: IPR value, sum(|psi|^4)
%
% Example:
%   ipr = calculate_ipr(eigenvector);

ipr = sum(abs(psi).^4);

end

