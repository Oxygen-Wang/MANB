function seq = fibonacci_numeric(n, S_A, S_B)
% Generate numeric array for Fibonacci sequence
%
% Parameters:
%   n: Number of substitutions
%   S_A: Numeric value corresponding to A
%   S_B: Numeric value corresponding to B
%
% Returns:
%   seq: Generated numeric sequence (column vector)
%
% Example:
%   seq = fibonacci_numeric(10, 1, 2);

char_seq = fibonacci_substitution(n);
seq = char_seq_to_numeric(char_seq, S_A, S_B);

end












