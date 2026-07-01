function overlap_ratio = calculate_overlap(interval1, interval2)
% Calculate overlap ratio between two intervals
%
% Parameters:
%   interval1: First interval [start, end]
%   interval2: Second interval [start, end]
%
% Returns:
%   overlap_ratio: Overlap ratio

start_overlap = max(interval1(1), interval2(1));
end_overlap = min(interval1(2), interval2(2));

if start_overlap >= end_overlap
    overlap_ratio = 0;
else
    overlap_length = end_overlap - start_overlap;
    min_length = min(interval1(2)-interval1(1), interval2(2)-interval2(1));
    overlap_ratio = overlap_length / min_length;
end

end













