function groups = contiguous_groups(idx_vec)
% Group contiguous indices
%
% Parameters:
%   idx_vec: Index vector
%
% Returns:
%   groups: Contiguous groups array [start, end] (N x 2)

if isempty(idx_vec)
    groups = [];
    return;
end
idx_vec = idx_vec(:)';
breaks = find(diff(idx_vec) > 1);
starts = [idx_vec(1), idx_vec(breaks+1)];
ends = [idx_vec(breaks), idx_vec(end)];
groups = [starts(:), ends(:)];

end













