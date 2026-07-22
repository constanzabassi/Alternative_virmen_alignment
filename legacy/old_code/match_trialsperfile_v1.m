function [file_estimated_trial_info,file_matching_trials,trial_times,file_digidata_trial_info] = match_trialsperfile(digidata_its, good_dataset, trial_info,sound_condition_array,alignment_info)
% find its in the data that best match the its for each trial dividing files into trials that match them
for file = 1:length(digidata_its)
file
end_trial_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps>.8*digidata_its(file).sync_sampling_rate));
start_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps>.8*digidata_its(file).sync_sampling_rate)+1);
% incorrect_trials_digidata_time = find(digidata_its(file).it_gaps >.8*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .95*digidata_its(file).sync_sampling_rate);
% correct_trials_digidata_time = find(digidata_its(file).it_gaps > .95*digidata_its(file).sync_sampling_rate);
start_trials_digidata_time =  digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate)+1);
end_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate));
big_gaps = find(digidata_its(file).it_gaps > .25*digidata_its(file).sync_sampling_rate); 
digidata_time.start_trials_digidata_time = start_trials_digidata_time;
digidata_time.end_trial_digidata_time = end_trial_digidata_time;
digidata_time.start_iti_digidata_time = start_iti_digidata_time; 
digidata_time.end_iti_digidata_time = end_iti_digidata_time;
trial_times(file) = digidata_time;
% figure(); hold on; plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(digidata_its(file).locs((big_gaps)),0,'*c'); plot(ex_data(:,7)); hold off
% distance between last correct big bag and next gap is  about 90 (could be
% as low as 74 probs high to 100

% distance between last correct big bag and next gap is  about 130 (ex
% values 159, 155, 140?, 132?

% for good datasets vs bad I am missing about half of the iterations so
% these difference between gaps has to be multiplied by 2 for good ones!
if good_dataset == 1
    good_dataset_conversion = 2;
else
    good_dataset_conversion = 1;
end

%write code to say if this big gap is followed by another big gap at this time point then!
difference_between_gaps = diff(big_gaps);
possible_outcomes = {};
count = 0;

for g = 2:length(difference_between_gaps)
    if difference_between_gaps(g) > 70*good_dataset_conversion && difference_between_gaps(g) < 101*good_dataset_conversion && digidata_its(file).it_gaps(big_gaps(g)) > .95*digidata_its(file).sync_sampling_rate
        count = count+1;
        possible_outcomes(count).correct = 1;
        [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g-1)+1) - alignment_info(file).frames_times));
        possible_outcomes(count).start_trial_digidata_time = digidata_its(file).locs(big_gaps(g-1)+1);
        if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
            possible_outcomes(count).frame_start_trial = closest_frame;
        else
            possible_outcomes(count).frame_start_trial = nan;
        end
        possible_outcomes(count).end_trial_digidata_time = digidata_its(file).locs(big_gaps(g)); %last iteration of current trial during maze (after this is the gap between end of maze and ITI)
       
        possible_outcomes(count).start_ITI_digidata_time = digidata_its(file).locs(big_gaps(g)+1); 
        
        [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g+1)) - alignment_info(file).frames_times));
        if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
            possible_outcomes(count).frame_end_ITI = closest_frame;
        else
            possible_outcomes(count).frame_end_ITI = nan;
        end
        possible_outcomes(count).end_ITI_digidata_time = digidata_its(file).locs(big_gaps(g+1));

    elseif difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion && digidata_its(file).it_gaps(big_gaps(g)) > .8*digidata_its(file).sync_sampling_rate && digidata_its(file).it_gaps(big_gaps(g)) < .95*digidata_its(file).sync_sampling_rate
        count = count+1;
        possible_outcomes(count).correct = 0;
        possible_outcomes(count).start_trial_digidata_time = digidata_its(file).locs(big_gaps(g-1)+1);
        [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g-1)+1) - alignment_info(file).frames_times));
        if val < .02*digidata_its(file).sync_sampling_rate
            possible_outcomes(count).frame_start_trial = closest_frame;
        else
            possible_outcomes(count).frame_start_trial = nan;
        end
        possible_outcomes(count).end_trial_digidata_time = digidata_its(file).locs(big_gaps(g));
        possible_outcomes(count).start_ITI_digidata_time = digidata_its(file).locs(big_gaps(g)+1); 

        possible_outcomes(count).end_ITI_digidata_time = digidata_its(file).locs(big_gaps(g+1));
        [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g+1)) - alignment_info(file).frames_times));
        if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
            possible_outcomes(count).frame_end_ITI = closest_frame;
        else
            possible_outcomes(count).frame_end_ITI = nan;
        end

        %TO DEAL WITH WEIRD TRIALS!!!!!! where sounds are played in the
        %wrong place or there are iterations missing from the end of the
        %trial
    elseif g>1 && min(abs(digidata_its(file).locs(big_gaps(g-1)+1) -digidata_time.start_trials_digidata_time)) == 0
            [val,trial_id] = min(abs(digidata_time.start_trials_digidata_time - digidata_its(file).locs(big_gaps(g-1)+1)));
            count = count+1;
            possible_outcomes(count).weird_trial = count;
            %assuming I am not loosing as many iterations during the ITI
            if difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion
                possible_outcomes(count).correct = 0;
            end

            if difference_between_gaps(g) > 70*good_dataset_conversion && difference_between_gaps(g) < 101*good_dataset_conversion
                possible_outcomes(count).correct = 1;
            end
            possible_outcomes(count).start_trial_digidata_time = digidata_time.start_trials_digidata_time(trial_id); 
            [val,closest_frame] = min(abs(digidata_time.start_trials_digidata_time(trial_id) - alignment_info(file).frames_times));
                if val < .02*digidata_its(file).sync_sampling_rate
                    possible_outcomes(count).frame_start_trial = closest_frame;
                else
                    possible_outcomes(count).frame_start_trial = nan;
                end
                [val,trial_id] = min(abs(digidata_time.end_trial_digidata_time - digidata_its(file).locs(big_gaps(g))));
                possible_outcomes(count).end_trial_digidata_time = digidata_time.end_trial_digidata_time(trial_id);
                [val,trial_id] = min(abs(digidata_time.start_iti_digidata_time - digidata_its(file).locs(big_gaps(g)+1)));
                possible_outcomes(count).start_ITI_digidata_time = digidata_time.start_iti_digidata_time(trial_id);

                [val,trial_id] = min(abs(end_iti_digidata_time - digidata_its(file).locs(big_gaps(g+1))));
                possible_outcomes(count).end_ITI_digidata_time = digidata_time.end_iti_digidata_time(trial_id);
                [val,closest_frame] = min(abs(digidata_time.end_iti_digidata_time(trial_id)- alignment_info(file).frames_times));
                if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
                    possible_outcomes(count).frame_end_ITI = closest_frame;
                else
                    possible_outcomes(count).frame_end_ITI = nan;
                end
%              elseif difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion && digidata_its(file).it_gaps(big_gaps(g)) > .95*digidata_its(file).sync_sampling_rate
%                 possible_incorrect_trial = [possible_incorrect_trial,digidata_its(file).locs(big_gaps(g))];
%                 count = count+1;
%                 possible_outcomes(count).correct = nan;
%                 possible_outcomes(count).start_trial_digidata_time = digidata_its(file).locs(big_gaps(g-1)+1);
%                 [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g-1)+1) - alignment_info(file).frames_times));
%                 if val < .02*digidata_its(file).sync_sampling_rate
%                     possible_outcomes(count).frame_start_trial = closest_frame;
%                 else
%                     possible_outcomes(count).frame_start_trial = nan;
%                 end
%                 possible_outcomes(count).end_trial_digidata_time = digidata_its(file).locs(big_gaps(g));
%                 possible_outcomes(count).start_ITI_digidata_time = digidata_its(file).locs(big_gaps(g)+1); 
%         
%                 possible_outcomes(count).end_ITI_digidata_time = digidata_its(file).locs(big_gaps(g+1));
%                 [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g+1)) - alignment_info(file).frames_times));
%                 if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
%                     possible_outcomes(count).frame_end_ITI = closest_frame;
%                 else
%                     possible_outcomes(count).frame_end_ITI = nan;
%                 end
    end
    
end
%check to make sure we are not skipping any trials!!
ex_start = start_trials_digidata_time;
ex_start(end) = []; %delete bc it would not be in the difference matrix
no_skips = ismember(ex_start,[possible_outcomes.start_trial_digidata_time]);
if sum(find(no_skips == 0))>1
    fprintf('Something is wrong, might be skipping trials\n');
else
    fprintf('Trials seem fine\n');
end

ex_data = abfload(strcat(digidata_its(file).directory));
figure(file);clf; hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));hold off; movegui(gcf,'center');

% combine sound condition and trial outcome based on digidata time!
% finds difference between start trial and sound offset
estimated_trial_info = [];weird_trials = []; count_t=0;trial_id = [];
for t = 1:length(possible_outcomes)
    %[val,closest_sound] = min(abs(possible_outcomes(t).digidata_time - sound_condition_array(file).file(:,3)));
    
    %below works but want to get rid of trials without frames at start and
    %end of session
%     if find(possible_outcomes(t).digidata_time - sound_condition_array(file).file(:,3)<0,1) > 1 %if it's 1 then the trial might be too short to look at
%         [closest_sound] = find(possible_outcomes(t).digidata_time - sound_condition_array(file).file(:,3)<0,1)-1; %find first value where it becomes negative and use the one before it
%         estimated_trial_info(t,:) = [possible_outcomes(t).correct;sound_condition_array(file).file(closest_sound,1)];
%         if sound_condition_array(file).file(closest_sound,3)> possible_outcomes(t).digidata_time  %if the sound plays during the ITI make it a NaN
%                 weird_trials = [weird_trials,t];
%         end
%     else
%         estimated_trial_info(t,:) = [possible_outcomes(t).correct;sound_condition_array(file).file(1,1)];
%     end

    if ~isnan(possible_outcomes(t).frame_start_trial)== 1%if it's 1 then the trial might be too short to look at
        count_t = count_t+1;
        [closest_sound] = find(possible_outcomes(t).end_trial_digidata_time - sound_condition_array(file).file(:,3)<0,1)-1; %find first value where it becomes negative and use the one before it
        estimated_trial_info(count_t,:) = [possible_outcomes(t).correct;sound_condition_array(file).file(closest_sound,1)];
        trial_id(count_t) = t;
        if sound_condition_array(file).file(closest_sound,2)< possible_outcomes(t).start_trial_digidata_time  %if the sound plays during the ITI make it a NaN
                weird_trials = [weird_trials,count_t];
                estimated_trial_info(count_t,:) = [possible_outcomes(t).correct;nan];
        end
    end
end
file_estimated_trial_info(file).estimated_trial_info = estimated_trial_info;
file_estimated_trial_info(file).trial_id = trial_id;
file_estimated_trial_info(file).weird_trials = weird_trials;
file_digidata_trial_info(file).estimated_digidata = possible_outcomes;
%figure(); hold on;plot(ex_data(:,6));plot(possible_correct_trial,ex_data(possible_correct_trial,6),'*c');plot(possible_incorrect_trial,ex_data(possible_incorrect_trial,6),'*r');hold off
%% estimate trials that most closely resemble each other

% Example data (replace with your actual data)
trueMatrix = [trial_info.correct; trial_info.condition]'; % Replace with your true trial information matrix
estimatedMatrix = estimated_trial_info; % Replace with your estimated trial information matrix

sectionSize = size(estimatedMatrix, 1);
numSections = size(trueMatrix, 1) - sectionSize + 1;
mseValues = zeros(numSections, 1);

for i = 1:numSections
    trueSection = trueMatrix(i:i+sectionSize-1, :);
    mseValues(i) = nanmean((trueSection - estimatedMatrix).^2, 'all');
end

[minMSE, bestTrueSectionIndex] = min(mseValues);

file_matching_trials(file,:) = [bestTrueSectionIndex, bestTrueSectionIndex + sectionSize - 1];
fprintf('Best matching true section of trials: %d to %d\n', bestTrueSectionIndex, bestTrueSectionIndex + sectionSize - 1);
fprintf('Minimum MSE value: %.4f\n', minMSE);
end
% verify that trials make sense... they go in chronological order
for file = 1:length(digidata_its) 
    if file > 1
        if file_matching_trials(file,1)>file_matching_trials(file-1,2) == 1
            fprintf('Trial order makes sense!\n'); %\n next line
        else
            fprintf('Error! Trials are out of order! Need to go back and check\n')
        end
    end
end