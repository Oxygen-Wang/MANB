function seq = thuemorse_substitution(n)
% Generate Thue-Morse substitution sequence (string form)
%
% Parameters:
%   n: Number of substitutions
%
% Returns:
%   seq: Generated character sequence (string)
%
% Description:
%   Thue-Morse substitution rules:
%   S0 = A
%   S_{n+1} = S_n + ~S_n (where ~ means swap A and B)
%
% Example:
%   seq = thuemorse_substitution(5);

S0 = 'A';

if n == 0
    seq = S0;
    return;
end

S_prev = S0;

for i = 1:n
    S_conj = strrep(S_prev, 'A', 'x');
    S_conj = strrep(S_conj, 'B', 'A');
    S_conj = strrep(S_conj, 'x', 'B');
    S_new = [S_prev S_conj];
    S_prev = S_new;
end

seq = S_new;

end












