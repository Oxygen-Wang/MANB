function seq = thuemorse_numeric(n, S_A, S_B)
% Generate numeric array for Thue-Morse sequence
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
%   seq = thuemorse_numeric(6, 1, 2);

char_seq = thuemorse_substitution(n);
seq = char_seq_to_numeric(char_seq, S_A, S_B);

end












