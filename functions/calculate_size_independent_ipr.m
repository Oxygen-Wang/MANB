function ipr = calculate_size_independent_ipr(psi)
% Calculate size-independent IPR (Inverse Participation Ratio)
%
% Parameters:
%   psi: Normalized eigenvector (row or column vector)
%
% Returns:
%   ipr: Size-independent IPR value, n * sum(|psi|^4)
%   where n is the system size (length of psi)
%
% Description:
%   Size-independent IPR is defined as: IPR = n * sum(|psi|^4)
%   where psi is normalized (norm(psi) = 1)
%
% Example:
%   psi = eigenvector / norm(eigenvector);
%   ipr = calculate_size_independent_ipr(psi);

% Ensure psi is a column vector
if isrow(psi)
    psi = psi';
end

% Get system size
n = length(psi);

% Calculate size-independent IPR
ipr = n * sum(abs(psi).^4);

end












