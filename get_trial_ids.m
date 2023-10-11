function [file_trial_ids,file_digidata_trial_info] = get_trial_ids(file_matching_trials,file_digidata_trial_info,file_estimated_trial_info)

for file = 1:size(file_matching_trials,1)
%% 1) find the first trial that would have imaging data and start with this trial if it is a full trial
frame_start_trial = [file_digidata_trial_info(file).estimated_digidata.frame_start_trial];
frame_end_ITI = [file_digidata_trial_info(file).estimated_digidata.frame_end_ITI];
trial_id = file_estimated_trial_info(file).trial_id; %within each file
trial_id(2,:) = [file_matching_trials(file,1):file_matching_trials(file,2)]; %relative to all other files
all_trials = 1:length(file_digidata_trial_info(file).estimated_digidata);

excluded_trials = sort(union(find(isnan(frame_start_trial)),find(isnan(frame_end_ITI)))); %no frames from start trial to end ITI within a trial
included_trials = setdiff(all_trials,excluded_trials);

if included_trials(1) == trial_id(1,1) %making sure I am not messing up the trial ids!
    start_trial = trial_id(2,find(trial_id(1,:) == included_trials(1)));
end

end_trial = trial_id(2,find(trial_id(1,:) == included_trials(end)));
file_trial_ids(file,:) = [start_trial,end_trial];

%put trials into context of all other trials
for t = all_trials
    if ismember(t,trial_id(1,:))
        file_digidata_trial_info(file).estimated_digidata(t).trial_id=trial_id(2,find(trial_id(1,:) ==t));
    else
        file_digidata_trial_info(file).estimated_digidata(t).trial_id = nan;
    end
end

end
%figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');

