function [rewards_per_file,reward_loc_pure_frames,reward_loc_end_trial,reward_loc_pure] = find_reward(virmen_it,sound_condition_array,mdl_end_trial_sol,mdl_pure_sol,file_trial_ids,trial_its,digidata_its,trial_info,alignment_info)
% for each trial find distance between end trial and start ITI in ms
rewards_per_file ={};
reward_loc_end_trial = zeros(length(trial_info),2);
reward_loc_pure = zeros(length(trial_info),4);
reward_loc_pure_frames = zeros(length(trial_info),4);
for file = 1:length(virmen_it) 
    %initialize variable to save
    it_distance = [];
    within_all_trials = file_trial_ids(file,1):file_trial_ids(file,2);
    count = 0;
    %loop across trials inside file and find distance between end of trial
    %and start of ITI
    
        %find frames based on pure tones
    min_distance = mode(diff(alignment_info(file).frame_times))+(2*alignment_info(file).sync_sampling_rate/1000); %minimum distance between frames in digidata time

    for trial = within_all_trials
        start_it = virmen_it(file).actual_it_values(1);
        end_it = trial_its.end_trial_its-start_it+1;
        start_iti_it = trial_its.start_iti_its-start_it+1;
        sound_trial = find([sound_condition_array(file).ITI_sounds(:,2)]>virmen_it(file).it_times(end_it(trial))& [sound_condition_array(file).ITI_sounds(:,2)]<virmen_it(file).it_times(start_iti_it (trial)));
        if trial_info(trial).correct == 1 %if mouse was rewarded
            count = count+1;
            temp = virmen_it(file).it_times(start_iti_it (trial)) - virmen_it(file).it_times(end_it(trial));
            temp = temp*1000/digidata_its(file).sync_sampling_rate; %convert to ms
            it_distance = [it_distance,temp];
            
            %use distance information to determine probable reward location
            %1) compare to distance from end of trial 
            reward_loc_end_trial(trial,1) = predict(mdl_end_trial_sol,it_distance(count)); %distance from end of trial iteration
            reward_loc_end_trial(trial,2) = virmen_it(file).it_times(end_it(trial)) + reward_loc_end_trial(trial,1)/1000*digidata_its(file).sync_sampling_rate; %distance in digidata time
            %2) compare to distance from pure tone 
            reward_loc_pure(trial,1) = predict(mdl_pure_sol,it_distance(count)); %distance from end of trial iteration
            reward_loc_pure(trial,2) = sound_condition_array(file).ITI_sounds(sound_trial,2) - reward_loc_pure(trial,1)/1000*digidata_its(file).sync_sampling_rate; %distance in digidata time
            reward_loc_pure(trial,3) = sound_condition_array(file).ITI_sounds(sound_trial,2);
            reward_loc_pure(trial,4) = sound_condition_array(file).ITI_sounds(sound_trial,3);
            %find frames
            reward_loc_pure_frames(trial,3) = find_closest_frames(alignment_info(file).frame_times,reward_loc_pure(trial,2),min_distance);
            reward_loc_end_trial(trial,3) = find_closest_frames(alignment_info(file).frame_times,reward_loc_end_trial(trial,2),min_distance);

            %saving frames of end of trial in this structure/ maybe more
            %accurate?
            reward_loc_pure_frames(trial,1) = find_closest_frames(alignment_info(file).frame_times,reward_loc_end_trial(trial,2),min_distance);
        end
        reward_loc_pure(trial,3) = sound_condition_array(file).ITI_sounds(sound_trial,2);
        reward_loc_pure(trial,4) = sound_condition_array(file).ITI_sounds(sound_trial,3);

        reward_loc_pure_frames(trial,2) = find_closest_frames(alignment_info(file).frame_times,reward_loc_pure(trial,3),min_distance);
        reward_loc_pure_frames(trial,3) = find_closest_frames(alignment_info(file).frame_times,reward_loc_pure(trial,4),min_distance);
        
    
    end
    %use distance information to determine probable reward location
    %1) compare to distance from end of trial

    rewards_per_file(file).it_distance = it_distance;

    
end
% figure plot the difference between the two predictions
difference =  reward_loc_pure(:,2)- reward_loc_end_trial(:,2);
difference = difference*1000/digidata_its(1).sync_sampling_rate;
figure(44);clf; plot(difference); xlabel('trial id'); ylabel('difference between pure and end trial prediction')
text(1,1,['mean difference: ' num2str(mean(abs(difference(difference~=0))))])
