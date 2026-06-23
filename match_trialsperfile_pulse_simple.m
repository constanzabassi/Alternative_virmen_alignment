function [file_digidata_trial_info,file_matching_trials] = match_trialsperfile_pulse_simple(digidata_its, trial_info,sound_condition_array,task_info,data)
%adjusted distances for iteration gaps! works for LR pulse maze

%initialize variable to save trials
file_matching_trials = [];

% find its in the data that best match the its for each trial dividing files into trials that match them
for file = 1:length(digidata_its)
file

%% get some digidata time information
%get indices without gaps to eliminate
gap_between_iti_start = 0.074; %.25
gap_betwee_end_iti = 0.1; %.8
end_trials_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps>gap_betwee_end_iti*digidata_its(file).sync_sampling_rate));
end_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps >gap_between_iti_start*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < gap_betwee_end_iti*digidata_its(file).sync_sampling_rate));

[~,indx] = intersect( digidata_its(file).locs,end_trials_digidata_time(:));
[~,indxiti] = intersect( digidata_its(file).locs,end_iti_digidata_time(:));

%use end ITI and end trial indices to get these
start_iti_digidata_time = digidata_its(file).locs(indx+1);
start_trials_digidata_time =  digidata_its(file).locs(indxiti+1);

file_digidata_trial_info(file).start_iti_digidata_time = start_iti_digidata_time;
file_digidata_trial_info(file).start_trials_digidata_time = start_trials_digidata_time;
file_digidata_trial_info(file).end_trials_digidata_time = end_trials_digidata_time;
file_digidata_trial_info(file).end_iti_digidata_time = end_iti_digidata_time;



ex_data = abfload(strcat(digidata_its(file).directory));
figure(55);clf;
title(strcat('Detected trial events file # ', num2str(file)));
hold on;aa = plot(ex_data(:,task_info.channel_number(1)));bb = plot(ex_data(:,task_info.channel_number(2)),'color',[0.7 0.7 0.7]);  cc = plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-b');
a = plot(end_trials_digidata_time,0,'*g');b = plot(start_trials_digidata_time,0,'*c');c= plot(end_iti_digidata_time,0,'*r'); movegui(gcf,'center');
%plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-r');
legend([aa bb cc  a(1) b(1) c(1)],'Imaging frames','Virmen its','Speaker 1', 'end trial', 'start trial', 'end iti')
if length(task_info.channel_number)>3
    dd = plot(rescale(ex_data(:,task_info.channel_number(4)),-1,0),'-m');
elseif length(task_info.channel_number)>4
    plot(rescale(ex_data(:,task_info.channel_number(5)),-1,0),'-r')
end
hold off;

% pause


%%
if ~isempty(digidata_its(file).pos_pks)%if positive peaks exist try to use them instead to assign trials to iterations

    [trial_its,trial_its_time] = virmen_it_rough_estimation(data); %get iterations that go with each trial
    assigned_pk_trial = find(digidata_its(file).pos_pks(1)*1e5 - trial_its.start_trial_its <= 0 ,1,'first'); 
    assgined_pk_time = digidata_its(file).pos_loc(1);
    within_file_trial = find(assgined_pk_time - file_digidata_trial_info(file).start_trials_digidata_time<= 0 ,1,'first'); 
    start_trial_index = assigned_pk_trial - within_file_trial;
    end_trial_index = start_trial_index + length(sound_condition_array(file).VR_sounds);%have to figure out if I wanna use start trial or end trial
    
    file_matching_trials(file,:) = [start_trial_index,end_trial_index];
    file_digidata_trial_info(file).trial_id = 1:length(end_trial_index - start_trial_index);


else
    
    % Example data (replace with your actual data)
    trueMatrix = [trial_info.correct; trial_info.condition]'; % Replace with your true trial information matrix
    estimatedMatrix = possible_outcome; % Replace with your estimated trial information matrix
    
    % Assuming trueMatrix and estimatedMatrix are appropriately defined
    sizeTrueMatrix = size(trueMatrix);
    trueConditions = cell(sizeTrueMatrix(1), 1);
    
    %Convert arrays into string - helps deal with multiple conditions
    for i = 1:sizeTrueMatrix(1)
        conditionStr = strjoin(cellstr(num2str(trueMatrix(i, 2))));
        correctness = trueMatrix(i, 1); % assuming 1st column contains correctness
        if isnan(correctness)
            correctnessStr = 'NaN';
        else
            correctnessStr = num2str(correctness);
        end
        trueConditions{i} = strcat(correctnessStr, '_', conditionStr);
    end
    
    size_estimatedMatrix = size(estimatedMatrix);
    estimatedConditions = cell(size_estimatedMatrix(1), 1);
    for i = 1:size_estimatedMatrix(1)
        conditionStr = strjoin(cellstr(num2str(estimatedMatrix{i, 2})));
        correctness = estimatedMatrix{i, 1};
        if isnan(correctness)
            correctnessStr = 'NaN';
        else
            correctnessStr = num2str(correctness);
        end
        estimatedConditions{i} = strcat(correctnessStr, '_', conditionStr);
    end
    
    %% match trials using levenshteinDistance (smallest change to strings gives distance)
    if file == 1
        startTrial = 1;
        bestTrialIndices = findBestMatchingTrials(trueConditions,estimatedConditions,startTrial);
    else
        startTrial = file_matching_trials(file-1,1);
        bestTrialIndices = findBestMatchingTrials(trueConditions,estimatedConditions,startTrial);
    end
    file_matching_trials(file,:) = [bestTrialIndices(1), bestTrialIndices(2)];
end
end
% verify that trials make sense... they go in chronological order
for file = 1:length(digidata_its) 
    if file > 1
        if file_matching_trials(file,1)>=file_matching_trials(file-1,2) == 1
            fprintf('Trial order makes sense!\n'); %\n next line
        else
            fprintf('Error! Trials are out of order! Need to go back and check\n')
            keyboard;
        end
    end
end