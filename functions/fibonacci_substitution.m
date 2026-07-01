function seq = fibonacci_substitution(n)
% Generate Fibonacci substitution sequence (string form)
%
% Parameters:
%   n: Number of substitutions
%
% Returns:
%   seq: Generated character sequence (string)
%
% Description:
%   Fibonacci substitution rules:
%   A -> AB
%   B -> A
%
% Example:
%   seq = fibonacci_substitution(5);

word = 'B';

for i = 1:n
    newWord = '';
    for j = 1:length(word)
        if word(j) == 'A'
            newWord = [newWord 'AB'];
        else
            newWord = [newWord 'A'];
        end
    end
    word = newWord;
end

seq = word;

end












