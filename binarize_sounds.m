function sounds_per_file = binarize_sounds(virmen_it,sound_condition_array, trial_its,sound_info,sound_st)
% sound_onsets_speakers = [sound_condition_array(file).VR_sounds{:,2}]; 
% sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
% sound_onsets_iterations = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > trial_its.start_trial_its(start_trial_number) & trial_its.sound_trigger_its <trial_its.end_iti_its(end_trial_number)))+6; %sound happens within 7 iterations
% possible_sound_onsets = virmen_it(file).it_times(sound_onsets_iterations);
overall_diff = [];
for file = 1:length(virmen_it)
    within_trials = virmen_it(file).file_trial_id_start:virmen_it(file).file_trial_id_end;
    within_trials_all = virmen_it(file).start_trial_number+1:virmen_it(file).end_trial_number;
    
    %sound_trigger_its_diff = virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1)-virmen_it(file).it_times(trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1-1);
    sound_trigger = trial_its.sound_trigger_its(within_trials_all)-virmen_it(file).actual_it_values(1)+1;
    sound_trigger_time = virmen_it(file).it_times(sound_trigger); %need to figure out indexing here
    
    sound_onsets = cell(1,1);
    sound_offsets= cell(1,1);
    [a,b] =find(sound_trigger_time(1) - [sound_condition_array(file).VR_sounds{:,3}]<0,1,'first');
    first_sound_trial = b-1;
    end_trial_time = virmen_it(file).it_times(trial_its.end_trial_its(within_trials_all-2+first_sound_trial)-virmen_it(file).actual_it_values(1)+1);
    first_onset = [];
    for s = 1:length(sound_trigger_time)
        if isnan(sound_condition_array(file).VR_sounds{s+first_sound_trial,3}) || ~isnan(sound_condition_array(file).VR_sounds{s+1,3})
            onset = [];
            %if nan determine onset and estimate where sounds might be
            %it before gap is about ~180ms from sound onset
            %on for 1 sec off for 200ms
            onset = sound_trigger_time(s)+virmen_it(file).mean_sound_distance; %conver to ms
            sound_onsets{s,:} = round(onset:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
            sound_offsets{s,:} = round(onset+sound_info.sound_duration:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
            if size(sound_onsets{s,:},2)>size(sound_offsets{s,:},2)% if last sound is unifinsehd
                sound_offsets{s,:} = [sound_offsets{s,:},round(end_trial_time(s))];
            end
            first_onset = [first_onset,sound_onsets{s,1}(1,1)];
        end
    end
    file
%     overall_diff = [overall_diff, [[sound_condition_array(file).VR_sounds{1:21,2}]-difference_test(1:21)]];
    
sounds_per_file(file).onsets = [sound_onsets{:,:}];%[[sound_st(file).file.onset],[sound_onsets{:,:}]];
sounds_per_file(file).offsets = [sound_offsets{:,:}];%[[sound_st(file).file.offset],[sound_offsets{:,:}]];
binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
for p = 1:length(sounds_per_file(file).offsets)
    binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1;
end

%visualize the binary sounds
ex_data = abfload(strcat(virmen_it(file).directory));
figure(555);clf;
title(strcat('Predicted binary sounds in file # ', num2str(file)));
hold on;cc = plot(rescale(ex_data(:,sound_info.spkr_channel_number(1)),0,1),'-b');dd = plot(rescale(ex_data(:,sound_info.spkr_channel_number(2)),0,1),'-m');
a = plot(binary_sounds,'-k');
%plot(rescale(ex_data(:,sound_info.spkr_channel_number(3)),-1,0),'-r');
legend([ cc dd  a(1) ],'Speaker 1','Speaker 2', 'binary sounds')
if length(sound_info.spkr_channel_number)>2
    plot(rescale(ex_data(:,sound_info.spkr_channel_number(3)),0,1),'-r')
end
hold off;
pause

%organize with actually detected sounds
sound_onsets = cell(1,1);
sound_offsets= cell(1,1);
for s = 1:length(sound_trigger_time)
    if isnan(sound_condition_array(file).VR_sounds{s+first_sound_trial,3}) 
        onset = [];
        %if nan determine onset and estimate where sounds might be
        %it before gap is about ~180ms from sound onset
        %on for 1 sec off for 200ms
        onset = sound_trigger_time(s)+virmen_it(file).mean_sound_distance; %conver to ms
        sound_onsets{s,:} = round(onset:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
        sound_offsets{s,:} = round(onset+sound_info.sound_duration:sound_info.distance_within_sounds+sound_info.sound_duration:end_trial_time(s));
        if size(sound_onsets{s,:},2)>size(sound_offsets{s,:},2)% if last sound is unifinsehd
            sound_offsets{s,:} = [sound_offsets{s,:},round(end_trial_time(s))];
        end
    end
end
sounds_per_file(file).onsets = [[sound_st(file).file.onset],[sound_onsets{:,:}]];
sounds_per_file(file).offsets = [[sound_st(file).file.offset],[sound_offsets{:,:}]];
binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
for p = 1:length(sounds_per_file(file).offsets)
    binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1;
end

sounds_per_file(file).binary_sounds = binary_sounds;

dist = [];
for s = 1:length(first_onset)
    temp = min(abs([sound_condition_array(file).VR_sounds{:,2}] - first_onset(s)));
    dist = [dist,temp];
end
sounds_per_file(file).distance_to_true = dist*1000/sound_info.sync_sampling_rate; %convert to ms

end
 mean_al = [];var_al =[];for i = 1:length(virmen_it); mean_al = [mean_al,mean(sounds_per_file(i).distance_to_true)];var_al = [var_al,var(sounds_per_file(i).distance_to_true)];end
 figure(556);clf; scatter(mean_al,var_al);xlabel('mean distance from true sound in ms'); ylabel('variance in  ms');


