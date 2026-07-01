function numeric_seq = char_seq_to_numeric(char_seq, S_A, S_B)
% Convert character sequence to numeric array
%
% Parameters:
%   char_seq: Character sequence (string)
%   S_A: Numeric value corresponding to character 'A'
%   S_B: Numeric value corresponding to character 'B'
%
% Returns:
%   numeric_seq: Numeric sequence (column vector)
%
% Example:
%   numeric_seq = char_seq_to_numeric('ABAB', 1, 2);

N = length(char_seq);
numeric_seq = zeros(N, 1);

for i = 1:N
    if char_seq(i) == 'A'
        numeric_seq(i) = S_A;
    elseif char_seq(i) == 'B'
        numeric_seq(i) = S_B;
    else
        error('Unknown character: %c', char_seq(i));
    end
end

end












