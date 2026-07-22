%function [gaps_to_eliminate]
% Define the minimum threshold in milliseconds
min_threshold = 3.8*digidata_its(file).sync_sampling_rate;

% List of gap differences (you provided this list)
gap_diffs = diff(digidata_its(file).locs(big_gaps));

% Initialize an array to store gaps to eliminate
gaps_to_eliminate = [];

% Iterate through gap differences to find gaps too close together
for i = 1:length(gap_diffs)-1
    if gap_diffs(i) < min_threshold
        % Check if the next gap is far from others
        if i == 1 || (i > 1 && gap_diffs(i-1) >= min_threshold)
            gaps_to_eliminate = [gaps_to_eliminate, i+1];
        else
            gaps_to_eliminate = [gaps_to_eliminate, i];
        end
    end
end

% Eliminate duplicate values from the gaps_to_eliminate array
gaps_to_eliminate = unique(gaps_to_eliminate);

% % Display the gaps to be eliminated
% disp('Gaps to be eliminated:');
% disp(gaps_to_eliminate);

%minimum threshold is ~3.8sec (distance between end trial and ITI is 4 sec)
% min_threshold = 3.8*1000;
% [bad_gaps] = find(diff(digidata_its(file).locs(big_gaps))<min_threshold);
% eliminate = [];
% diff_bad_gaps = diff(bad_gaps);
% for g = 1:length(bad_gaps)-1 
%     if diff_bad_gaps(g)>1
%         eliminate = [eliminate,bad_gaps(g)+1];
%     end
% end
%want to eliminate #12, 34, 39