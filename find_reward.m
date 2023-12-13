function rewards_per_file = find_reward(virmen_it,sound_condition_array,mdl_end_trial_sol,mdl_pure_sol,file_trial_ids,trial_its)
% for each trial find distance between end trial and start ITI in ms
rewards_per_file ={};
for file = 1:length(virmen_it) 
    %initialize variable to save
    it_distance = [];
    %loop across trials inside file
    for trial = file_trial_ids(file,1):file_trial_ids(file,2)
        start_it = virmen_it(file).actual_it_values(1);
        end_it = trial_its.end_trial_its-start_it+1;
        start_iti_it = trial_its.start_iti_its-start_it+1;
        temp = virmen_it(file).it_times(start_iti_it (trial)) - virmen_it(file).it_times(end_it(trial));
        it_distance = [it_distance,temp];
    end
    rewards_per_file(file).it_distance = it_distance;
end
