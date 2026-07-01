function seq = doubleperiod_substitution(n)
% Generate Double Period substitution sequence (string form)
%
% Parameters:
%   n: Number of substitutions
%
% Returns:
%   seq: Generated character sequence (string)
%
% Description:
%   Double Period substitution rules:
%   S0 = A
%   S1 = AB
%   A -> AB
%   B -> AA
%
% Example:
%   seq = doubleperiod_substitution(5);

S0 = 'A';
S1 = 'AB';

if n == 0
    seq = S0;
    return;
elseif n == 1
    seq = S1;
    return;
end

S_prev = S1;

for i = 2:n
    newWord = '';
    for j = 1:length(S_prev)
        if S_prev(j) == 'A'
            newWord = [newWord 'AB'];
        elseif S_prev(j) == 'B'
            newWord = [newWord 'AA'];
        end
    end
    S_prev = newWord;
end

seq = S_prev;

end












