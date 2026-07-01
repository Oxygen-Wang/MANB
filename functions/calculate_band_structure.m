function [real_eigval_pro, imag_eigval_pro, eigvec_pro1, eigvec_pro2] = ...
    calculate_band_structure(k_pro, n_pro, dz, rho, c)
% Calculate band structure (eigenvalues and eigenvectors) for acoustic waveguide
%
% Parameters:
%   k_pro: k sampling point array (1 x num_k row vector)
%   n_pro: Waveguide cross-section sequence (num_n x 1 column vector)
%   dz: Length interval
%   rho: Density (kg/m^3)
%   c: Sound speed (m/s)
%
% Returns:
%   real_eigval_pro: Real part of eigenvalues (2 x num_k), two bands
%   imag_eigval_pro: Imaginary part of eigenvalues (2 x num_k), two bands
%   eigvec_pro1: Eigenvectors of first band (2 x num_k)
%   eigvec_pro2: Eigenvectors of second band (2 x num_k)
%
% Description:
%   For each k value, calculate transfer matrix, then solve for eigenvalues and eigenvectors.
%   Eigenvalues are processed and sorted to ensure correct band correspondence.
%
% Example:
%   [real_eigval, imag_eigval, eigvec1, eigvec2] = ...
%       calculate_band_structure(k_pro, n_pro, dz, rho, c);

num_k = length(k_pro);
num_n = length(n_pro);

% Initialize storage arrays
real_eigval_pro = zeros(2, num_k);
imag_eigval_pro = zeros(2, num_k);
eigvec_pro1 = zeros(2, num_k);  % Store eigenvectors of first band
eigvec_pro2 = zeros(2, num_k);  % Store eigenvectors of second band

% Main loop: calculate transfer matrix and solve for eigenvalues and eigenvectors
for ii = 1:num_k
    k = k_pro(ii);
    
    % Calculate transfer matrix
    TM_T = eye(2);
    for jj = 1:num_n
        S_curr = n_pro(jj);
        T1 = acoustic_TM(S_curr, k, dz, rho, c);
        TM_T = T1 * TM_T;
    end
    
    % Calculate eigenvalues and eigenvectors
    [eigvec0, eigval0] = eig(TM_T);
    val1 = -1i*log(eigval0(1,1));
    val2 = -1i*log(eigval0(2,2));
    
    % Classify bands according to eigenvalue characteristics
    if (abs(real(val1)-pi) < 1e-5) || (real(val1) < 1e-5)
        if imag(val1) > imag(val2)
            eigvec_pro1(:,ii) = eigvec0(:,1);
            eigvec_pro2(:,ii) = eigvec0(:,2);
        else
            eigvec_pro1(:,ii) = eigvec0(:,2);
            eigvec_pro2(:,ii) = eigvec0(:,1);
        end
    else
        if real(val1) > real(val2)
            eigvec_pro1(:,ii) = eigvec0(:,1);
            eigvec_pro2(:,ii) = eigvec0(:,2);
        else
            eigvec_pro1(:,ii) = eigvec0(:,2);
            eigvec_pro2(:,ii) = eigvec0(:,1);
        end
    end
    
    % Calculate and process real and imaginary parts of eigenvalues
    real_eigvals = diag(real(-1i*log(eigval0)));
    imag_eigvals = diag(imag(-1i*log(eigval0(1:2,1:2))));
    imag_copy = imag_eigvals;
    
    % Map real part of eigenvalues to [-pi, pi] interval
    real_eigvals = mod(real_eigvals + pi, 2*pi) - pi;
    real_copy = real_eigvals;
    
    % Sort real part of eigenvalues
    if real_eigvals(2) > real_eigvals(1)
        real_eigvals(1) = real_copy(2);
        real_eigvals(2) = real_copy(1);
    end
    real_eigvals(1) = -real_eigvals(2);
    
    % Sort imaginary part of eigenvalues
    if imag_eigvals(2) > imag_eigvals(1)
        imag_eigvals(1) = imag_copy(2);
        imag_eigvals(2) = imag_copy(1);
    end
    
    % Store results
    real_eigval_pro(:,ii) = real_eigvals;
    imag_eigval_pro(:,ii) = imag_eigvals;
end

end












