% [sound_outputs,trialConditions, condition_onset_array_all]=
function [sound_outputs,trialConditions, condition_onset_array_all]= find_spkr_output_task_new(server,mouse,date,alignment_info,spkr_channel_number,string,detection_threshold,distance_between_sounds,distance_within_sounds,sound_duration) 
%pc=1 if windows, any other number if mac
cd(strcat(server,'\Connie\RawData\',num2str(mouse),'\wavesurfer\',num2str(date)));
sync_dir = dir(strcat('*',string,'*.abf'));
num_files = length(sync_dir);
file_ind = 0;

for file = 1:num_files
    file_ind = file_ind +1
    [sync_data,sync_sampling_interval,~]  = abfload(sync_dir(file).name);
    sync_sampling_rate = 1/sync_sampling_interval*1e6;
    sync_data=sync_data';
    sync_data=double(sync_data);
    rawSounds=sync_data(spkr_channel_number,:);
    frames_times{file} = alignment_info(file).frame_times; 
    rescaled_sounds=[];
    for s = 1:size(rawSounds,1)
        rescaled_sounds(s,:) = rescale(rawSounds(s,:),-1,1);
    end


    detection_threshold = 0.45;
    bin_sound_signal=[];
    pure_tone_signal=[];
    for s = 1:size(rawSounds,1)
    [binary_sound_signal,pure_tone] = process_sound_signal(rescaled_sounds(s,:),detection_threshold,sync_sampling_rate,15); %last value is smoothing factor
    bin_sound_signal(s,:) = binary_sound_signal;
    pure_tone_signal(s,:) = pure_tone;
    end

    originalVector = [];
    originalVector = bin_sound_signal;
    % Reverse 0s and 1s
    reversedSoundVector = 1 - originalVector;
    % Include all sounds in one vector
    soundVector = max(reversedSoundVector);

    %figure(); hold on;plot(soundVector); plot(rescaled_sounds(1,:));plot(rescaled_sounds(2,:)); hold off
    % Initialize a counter
    counter = 0;
    
    % Create a vector with numbers indicating the order of sound
    orderedVector = zeros(size(originalVector));
    
    groupedVector = [];
    % Iterate through the original vector
    for i = 1:length(soundVector)
        if soundVector(i) == 1
            if i == 1 || soundVector(i - 1) == 0 && soundVector(i - 2) == 0% -2 to deal with signal changing right after the offset of a sound to make sure it is included within the first group
                % Start of a new group
                counter = counter + 1;
            end
            groupedVector(i) = counter;
        end
    end
    %force the signal to go to zero at the end just to get the offset of
    %the final sound
    groupedVector(1,end-1:end) = 0;

    sounds1 = groupedVector;
    diff_sounds = diff(sounds1);
    [pks,locs] = findpeaks(abs(diff_sounds));
    sound_onset = []; sound_offset = [];
    sound_onset = locs(find(diff_sounds(locs)>0));
    sound_offset = locs(find(diff_sounds(locs)<0));
        % Print information for debugging
    fprintf('Total detected onsets: %d\n', length(sound_onset));
    fprintf('Total detected offsets: %d\n', length(sound_offset));
     % Find unpaired onsets and offsets
    unpaired_onset = setdiff(sounds1(sound_onset + 1), sounds1(sound_offset)); % onset not in offset
    unpaired_offset = setdiff(sounds1(sound_offset), sounds1(sound_onset + 1)); % offset not in onset
    
    onset_index = setdiff(sounds1(sound_onset+1),unpaired_onset);
    offset_index = setdiff(sounds1(sound_offset),unpaired_offset);
    ss = 0; sound_pairs = [];true_sound_pairs = [];
    for s = 1:min([length(onset_index),length(offset_index)])%min([length(sound_onset),length(sound_offset)])
        onset_sound = sound_onset(sounds1(sound_onset+1) == onset_index(s));
        offset_sound = sound_offset(sounds1(sound_offset) == offset_index(s));

        %adding this to deal with small changes in the signal right at onset or offset to make sure
        %I am using correct onset and offset
        if length(offset_sound)>1 && offset_sound(1) - onset_sound <10 %if the distance is large then first one is most likely the correct one
            sound_pairs(s,:) = [onset_sound;offset_sound(2)];
        else
            sound_pairs(s,:) = [onset_sound;offset_sound(1)];
        end
            if sound_pairs(s,1) < size(rawSounds,2) && sound_pairs(s,2) < size(rawSounds,2) %unsure why this would be bigger than array but does happen
                ss = ss+1;
                true_sound_pairs(ss,:) = sound_pairs(s,:);
            else
                sound_pairs(s,:) = nan;
            end
    end

    % Print the valid pairs for debugging
    fprintf('Number of valid pairs: %d\n', length(true_sound_pairs));
    disp(true_sound_pairs);
    
    %distance_between = 0.22; %distance between sounds within a trial
        
    difference = [true_sound_pairs(:,2) - true_sound_pairs(:,1)]; 
    all_trial_sounds = [];
    all_trial_sounds = true_sound_pairs(find(difference >sound_duration(1) & difference < sound_duration(2)),:); %sounds that are outside limits of sound duration
    % adding code to also include sounds that are cut off early
    count = 0; unfinished_sounds = [];unfinished_sounds_toadd =[];
    
    unfinished_sounds = setdiff(1:length(difference),find(difference >sound_duration(1) & difference < sound_duration(2)));
    if ~isempty(unfinished_sounds)
        for es = 1:length(unfinished_sounds)
            extra_sound = unfinished_sounds(es);
            if extra_sound > 1 && extra_sound<unfinished_sounds(end) && [sound_pairs(extra_sound,1) - sound_pairs(extra_sound-1,2)] < (distance_within_sounds*.1+distance_within_sounds) && [sound_pairs(extra_sound,1) - sound_pairs(extra_sound-1,2)] > (distance_within_sounds-(distance_within_sounds*.1))
                count = count+1;
                    unfinished_sounds_toadd(count,:) = [sound_pairs(extra_sound,:)];
            end
        end
        all_trial_sounds = sort([all_trial_sounds; unfinished_sounds_toadd]);
    end

    %classify sounds
     [sound_struc, condition_array, onset_array, offset_array,classified_sounds] = classify_sound_2spkr (reversedSoundVector,all_trial_sounds);

     numSounds = length(onset_array);
numConditions = length(unique(condition_array));
expectedDistance = distance_within_sounds; %distance_between*sync_sampling_rate; 

% Initialize cell arrays to store trial information for each condition
trialsPerCondition = cell(numConditions, 1);
groupNum = 0; % Initialize group number
onsettimes= [];
condition_group_array = [];
% Iterate over each sound
for i = 1:numSounds
    onsetTime = onset_array(i);
    offsetTime = offset_array(i);
    condition = condition_array(i);
    % Initialize or get the current condition's trials
    if isempty(trialsPerCondition{condition})
        trialsPerCondition{condition} = [onsetTime, offsetTime];
    else
        lastTrial = trialsPerCondition{condition}(end, :);
        lastOffsetTime = lastTrial(2);
        
        % Check if the current sound can be added to the last trial
        if abs(lastOffsetTime - onsetTime) <= expectedDistance+expectedDistance*.1
            % Extend the last trial
            trialsPerCondition{condition}(end, 2) = offsetTime;
        else
            % Create a new trial
            trialsPerCondition{condition} = [trialsPerCondition{condition}; [onsetTime, offsetTime]];
        end
    end
    % add to sound struct
    if i > 1 && onsetTime - offset_array(i - 1) <= distance_between_sounds
        % Assign the same group number as the previous sound
        sound_struc(i).trial_num = sound_struc(i - 1).trial_num;
        
    else
        % Increment the group number and assign to the current sound
        groupNum = groupNum + 1;
        sound_struc(i).trial_num = groupNum;
        onsettimes = [onsettimes,onsetTime];
        condition_group_array(groupNum, 2) = onsettimes(1);
    end
    
    % Store condition and group number information
    condition_group_array(groupNum, 1) = condition;
    
    condition_group_array(groupNum, 3) = offsetTime;
    condition_group_array(groupNum, 4) = groupNum;
    onsettimes=[];
end



%% Display the resulting trials

figure(110);clf
subplot(2,1,1)
% Plot the sound_array as a binary plot
plot(classified_sounds, 'b');
ylim([-0.5, 2.5]);
yticks([0, 1]);
yticklabels({'No Sound', 'Sound'});
xlabel('Time');
title('Sound Array');

hold on;

% Plot the trials for each condition
for c = 1:numConditions
    currentConditionTrials = trialsPerCondition{c};
    if ~isempty(currentConditionTrials)
        for t = 1:size(currentConditionTrials, 1)
            trial = currentConditionTrials(t, :);
            rectangle('Position', [trial(1), c-0.4, trial(2)-trial(1), 0.8], 'FaceColor', [1, 0, 0, 0.3], 'EdgeColor', 'none');
        end
    end
end

% Adjust y-axis and labels for the conditions
yticks(1:numConditions);
con_labels = {};for labels = 1:numConditions; con_labels(labels,:) = {strcat('Condition ',num2str(labels))};end
yticklabels(con_labels);%{'Condition 1', 'Condition 2'}); % Add labels as needed

hold off;
subplot(2,1,2)
colors_s = [0.7 0.3 0.9;0.4 0.6 0.9;0.9 0 0.6; 1 0.6 2];
hold on
% Plot the trials for each condition
for c = 1:numConditions
    currentConditionTrials = trialsPerCondition{c};
    if ~isempty(currentConditionTrials)
        for t = 1:size(currentConditionTrials, 1)
            trial = currentConditionTrials(t, :);
            rectangle('Position', [trial(1), 1-0.4, trial(2)-trial(1), 0.8], 'FaceColor', [1, 0, 0, 0.3], 'EdgeColor', 'none');
        end
    end
end
for s = 1:size(rescaled_sounds,1)
    plot(rescaled_sounds(s,:),'color',colors_s(s,:));
end
    plot(soundVector,'b')
hold off
xlabel('Normalized sounds');
title('Sound Array');
legend(con_labels)
%% outputs: sound struc has all sound onset and offsets for each individual sound
sound_outputs(file).file = sound_struc;
trialConditions(file).file = trialsPerCondition;
condition_onset_array_all(file).file = condition_group_array;
pause

end
