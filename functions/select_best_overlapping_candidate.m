function best_candidate = select_best_overlapping_candidate(candidates, reference_candidate)
% Select candidate with maximum slope that overlaps with reference candidate
%
% Parameters:
%   candidates: Candidate narrow band array
%   reference_candidate: Reference candidate [k_min, k_max, ...]
%
% Returns:
%   best_candidate: Best candidate narrow band

if isempty(candidates)
    best_candidate = [];
    return;
end

% Calculate overlap between each candidate and reference candidate
overlap_scores = zeros(size(candidates,1), 1);
for i = 1:size(candidates,1)
    overlap_scores(i) = calculate_overlap(reference_candidate(1:2), candidates(i,1:2));
end

% Only consider candidates with overlap
valid_mask = overlap_scores > 0;
if ~any(valid_mask)
    best_candidate = [];
    return;
end

valid_candidates = candidates(valid_mask, :);
valid_scores = overlap_scores(valid_mask);

% Select candidate with maximum average slope
[~, best_idx] = max(valid_candidates(:,4));
best_candidate = valid_candidates(best_idx, :);

end













