function IPR = calculate_ipr_for_all_states(eigenvectors)
% Calculate IPR values for all eigenstates
%
% Parameters:
%   eigenvectors: Eigenvector matrix, each column is an eigenvector (N x N matrix)
%
% Returns:
%   IPR: Array of IPR values for each eigenstate (N x 1 column vector)
%
% Example:
%   IPR = calculate_ipr_for_all_states(V);

N = size(eigenvectors, 2);
IPR = zeros(N, 1);

for k = 1:N
    psi = eigenvectors(:, k);
    IPR(k) = calculate_ipr(psi);
end

end

