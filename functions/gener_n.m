function n_t = gener_n(num_n, N, Sa, Sb)
% Generate periodic sequence
%
% Parameters:
%   num_n: Total number of elements
%   N: Number of periods
%   Sa: Value for first half of period
%   Sb: Value for second half of period
%
% Returns:
%   n_t: Generated sequence (num_n x 1 column vector)
%
% Example:
%   n_t = gener_n(100, 5, 0.01, 0.02);

n_t = zeros(num_n, 1);

period = num_n / N;

% Check if period is an integer
if mod(period, 1) ~= 0
    error('num_n / N must be an integer');
end

shift = round(period/4);

for ii = 0:N-1
    for jj = 0:period-1
        idx = mod(round(ii*period + jj + shift), num_n) + 1;
        if jj < period/2
            n_t(idx) = Sa;
        else
            n_t(idx) = Sb;
        end
    end
end

end

