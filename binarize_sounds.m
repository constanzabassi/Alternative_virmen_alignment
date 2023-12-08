function sounds_per_file = binarize_sounds(virmen_it,sound_condition_array, trial_its,sound_info,sound_st)
% sound_onsets_speakers = [sound_condition_array(file).VR_sounds{:,2}]; 
% sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
% sound_onsets_iterations = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > trial_its.start_trial_its(start_trial_number) & trial_its.sound_trigger_its <trial_its.end_iti_its(end_trial_number)))+6; %sound happens within 7 iterations
% possible_sound_onsets = virmen_it(file).it_times(sound_onsets_iterations);
overall_diff = [];
for file = 1:length(virmen_it)
    within_trials = virmen_it(file).file_trial_id_start:virmen_it(file).file_trial_id_start;
    within_trials_all = virmen_it(file).start_trial_number+1:virmen_it(file).end_trial_number;
    
    %sound_trigger_its_diff = virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1)-virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1-1);
    sound_trigger = trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1;
    end_trial_time = virmen_it(file).it_times(trial_its.end_trial_its(within_trials_all-1)-virmen_it(file).actual_it_values(1)+1);
    sound_trigger_time = virmen_it(file).it_times(sound_trigger); %need to figure out indexing here
    
    sound_onsets = cell(1,1);
    sound_offsets= cell(1,1);
    for s = 1:length(sound_trigger_time)
        if isnan(sound_condition_array(file).VR_sounds{s,3})
            onset = [];
            %if nan determine onset and estimate where sounds might be
            %it before gap is about ~180ms from sound onset
            %on for 1 sec off for 200ms
            onset = sound_trigger_time(s)+96.5*(1000/sound_info.sync_sampling_rate); %conver to ms
            sound_onsets{s,:} = onset:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s);
            sound_offsets{s,:} = onset+sound_info.sound_duration:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s);
            if size(sound_onsets{s,:},2)>size(sound_offsets{s,:},2)% if last sound is unifinsehd
                sound_offsets{s,:} = [sound_offsets{s,:},end_trial_time(s)];
            end
        end
    end
    file
    missing_sounds(file).onset = sound_onsets;
    missing_sounds(file).offset = sound_offsets;
    difference_test = sound_trigger_time+96.5*(1000/sound_info.sync_sampling_rate);
%     overall_diff = [overall_diff, [[sound_condition_array(file).VR_sounds{1:21,2}]-difference_test(1:21)]];
    
    sounds_per_file(file).onsets = [[sound_st(file).file.onset],sound_onsets{:,1}];
    sounds_per_file(file).offsets = [[sound_st(file).file.offset],sound_offsets{:,1}];
    binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
    for p = 1:length(sounds_per_file(file).offsets)
        binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1;
    end
    sounds_per_file(file).binary_sounds = binary_sounds;

end
figure();
plot(overall_diff);

