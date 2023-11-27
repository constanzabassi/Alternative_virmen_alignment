function [file_digidata_trial_info,file_matching_trials] = match_trialsperfile(digidata_its, good_dataset, trial_info,sound_condition_array,alignment_info)
%initialize variable to save trials
file_matching_trials = [];

% find its in the data that best match the its for each trial dividing files into trials that match them
for file = 1:length(digidata_its)
file
big_gaps = find(digidata_its(file).it_gaps > .25*digidata_its(file).sync_sampling_rate); 

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

% get rid of gaps too close together
% Define the minimum threshold in milliseconds
min_threshold = 3.8*digidata_its(file).sync_sampling_rate; %distance would be 4 sec
max_threshold = 6.0 * digidata_its(file).sync_sampling_rate; %max distance would be 6 sec
%%
% List of gap differences 
gap_diffs = diff(digidata_its(file).locs(big_gaps));

%Use sound information to see if gaps are about where they should be 
%keep in mind that if it's 8 locs I might be missing some sounds
%iterate through gaps
probable_end_iti_gaps = [];
probable_end_trial_gaps=[];
gaps_to_eliminate = [];
good_g = [];
for g = 1:length(digidata_its(file).locs(big_gaps))
%1) use ITI sounds to find gaps that are near the start of the next trial
difference = digidata_its(file).locs(big_gaps(g)) - sound_condition_array(file).ITI_sounds(:,3);
[a,~]=find(difference>0,1,'last'); %finds closest gap
%temp2 = [temp2,difference(a)]; % actual difference
% seems like if it's correct the difference is between 3.1s and 3.2s
% if it's incorrect the differencec is between 5.2s and maybe 5.3s?
if  any(sound_condition_array(file).ITI_sounds(a,1)==1) && any(difference(a) > 3*digidata_its(file).sync_sampling_rate) && any(difference(a) < 3.3*digidata_its(file).sync_sampling_rate)
    probable_end_iti_gaps(a,1) = [digidata_its(file).locs(big_gaps(g))];
    probable_end_iti_gaps(a,2) = [sound_condition_array(file).ITI_sounds(a,1)];
    good_g = [good_g,g];
elseif  any(sound_condition_array(file).ITI_sounds(a,1)==0 ) && any(difference(a) > 5*digidata_its(file).sync_sampling_rate) && any(difference(a) < 5.5*digidata_its(file).sync_sampling_rate)
    probable_end_iti_gaps (a,1) =[digidata_its(file).locs(big_gaps(g))];
    probable_end_iti_gaps(a,2) = [sound_condition_array(file).ITI_sounds(a,1)];%probable_end_iti_gaps(a+1,2) = [sound_condition_array(file).ITI_sounds(a+1,1)];
    good_g = [good_g,g];
end
%2) find sounds that are close to the end of the trial and compare to gaps
%that should also be close to the end of the trial
difference_VR = abs(digidata_its(file).locs(big_gaps(g)) - [sound_condition_array(file).VR_sounds{:,3}]);
difference_VR_noabs = (digidata_its(file).locs(big_gaps(g)) - [sound_condition_array(file).VR_sounds{:,3}]);
[~,b]=min((difference_VR)); %finds closest gap
if  any(difference_VR(b) < .2*digidata_its(file).sync_sampling_rate) %might need to change to 0.2 since that might be the greatest difference
%     if difference_VR_noabs(b) < 0 %use sound to say this is where the last gap happens
%         probable_end_trial_gaps(b,1) = sound_condition_array(file).VR_sounds(b,3);
%     else
%         probable_end_trial_gaps(b,1) = [digidata_its(file).locs(big_gaps(g))];
%     end

    probable_end_trial_gaps(b,1) = [digidata_its(file).locs(big_gaps(g))];
    probable_end_trial_gaps(b,2) = [sound_condition_array(file).VR_sounds{b,1}];
    good_g = [good_g,g];
end

%find difference between gaps to get rid of gaps that are too close
%together and don't meet the criteria above
if  ismember(g,good_g) && g<length(gap_diffs) && gap_diffs(g) < min_threshold
    if good_g(end)<length(digidata_its(file).locs(big_gaps)) 
    gaps_to_eliminate = [gaps_to_eliminate,digidata_its(file).locs(big_gaps(g+1))];
    end
end
end


%% get some digidata time information

%get indices without gaps to eliminate
end_trials_digidata_time = setdiff(digidata_its(file).locs(find(digidata_its(file).it_gaps>.8*digidata_its(file).sync_sampling_rate)),gaps_to_eliminate);
end_iti_digidata_time = setdiff(digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate)),gaps_to_eliminate);%digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate));

unpaired = {};
unpaired.iti = setdiff(end_iti_digidata_time,probable_end_iti_gaps(:,1));
unpaired.trial = setdiff(end_trials_digidata_time,probable_end_trial_gaps(:,1));
% pair up the probable end trial and end iti gaps
% pair up so end trial goes first and then end iti
trial_pairs = zeros(min(length(probable_end_trial_gaps),length(probable_end_iti_gaps)),2);
possible_outcome = [];
for t = 1:min(length(probable_end_trial_gaps),length(probable_end_iti_gaps))
    if probable_end_trial_gaps(t) < probable_end_iti_gaps(t) && probable_end_trial_gaps(t) > 0 %%assumes that you always start with the end trial which might not always be the case!
        trial_pairs(t,:) = [probable_end_trial_gaps(t),probable_end_iti_gaps(t)];
    elseif probable_end_trial_gaps(t) >0
        difference = end_iti_digidata_time-probable_end_trial_gaps(t); 
        [~,a]=find(difference>5.750*digidata_its(file).sync_sampling_rate & difference < 5.97*digidata_its(file).sync_sampling_rate); %finds closest gap incorrect will be ~5.8sec and correct should be ~4sec
        [~,b]=find(difference>3.85*digidata_its(file).sync_sampling_rate & difference < 4.05*digidata_its(file).sync_sampling_rate); %finds closest gap incorrect will be ~5.8sec and correct should be ~4sec
        trial_pairs(t,:) = [probable_end_trial_gaps(t), end_iti_digidata_time(a),end_iti_digidata_time(b)];
    elseif probable_end_iti_gaps(t) >0
        difference = probable_end_iti_gaps(t) - end_trials_digidata_time; 
        [~,a]=find(difference>5.750*digidata_its(file).sync_sampling_rate & difference < 5.97*digidata_its(file).sync_sampling_rate); %finds closest gap incorrect will be ~5.8sec and correct should be ~4sec
        [~,b]=find(difference>3.85*digidata_its(file).sync_sampling_rate & difference < 4.05*digidata_its(file).sync_sampling_rate); %finds closest gap incorrect will be ~5.8sec and correct should be ~4sec
        trial_pairs(t,:) = [end_trials_digidata_time(a),end_trials_digidata_time(b),probable_end_iti_gaps(t)];
    elseif t>1 && t<min(length(probable_end_trial_gaps),length(probable_end_iti_gaps)) % if both probable are zero then there is no sound nearby and need to figure it out on it's own OR call it a nan?
       difference = unpaired.trial - trial_pairs(t-1,1);
       [~,a]=find(difference>8*digidata_its(file).sync_sampling_rate ); 
       difference2 = unpaired.iti - trial_pairs(t-1,2);
       [~,b]=find(difference2>8*digidata_its(file).sync_sampling_rate ); 
       trial_pairs(t,:) = [unpaired.trial(a(1)),unpaired.iti(b(1))];
    else
        trial_pairs(t,:) =[nan,nan];
    end
end

end_trials_digidata_time = trial_pairs(:,1);
[~,indx] = intersect( digidata_its(file).locs,end_trials_digidata_time(:));
end_iti_digidata_time = trial_pairs(:,2);
[~,indxiti] = intersect( digidata_its(file).locs,end_iti_digidata_time(:));

%use end ITI and end trial indices to get these
start_iti_digidata_time = digidata_its(file).locs(indx+1);
start_trials_digidata_time =  digidata_its(file).locs(indxiti+1);

file_digidata_trial_info(file).start_iti_digidata_time = start_iti_digidata_time;
file_digidata_trial_info(file).start_trials_digidata_time = start_trials_digidata_time;
file_digidata_trial_info(file).end_trials_digidata_time = end_trials_digidata_time;
file_digidata_trial_info(file).end_iti_digidata_time = end_iti_digidata_time;

% build outcome array based on sounds and ITI sounds!
possible_outcome = {};
count_t = 0;trial_id = []; weird_trials = [];
for t = 1:length(sound_condition_array(file).VR_sounds)
     %test for weird trials where sound plays during ITI (virmen bug)
    if find(sound_condition_array(file).VR_sounds{t,2} > file_digidata_trial_info(file).start_iti_digidata_time & sound_condition_array(file).VR_sounds{t,2} <file_digidata_trial_info(file).end_iti_digidata_time & (file_digidata_trial_info(file).end_iti_digidata_time - file_digidata_trial_info(file).start_iti_digidata_time)<5.750*digidata_its(file).sync_sampling_rate)
        weird_trials = [weird_trials,t];
        file_digidata_trial_info(file).weird_trial = weird_trials;
    end
    %could change possible outcome for weird trials if the matching doesnt
    %work- seems to still work even with these incorrectly labeled trials
    possible_outcome{t,1} = [sound_condition_array(file).ITI_sounds(t,1)]; % correct vs incorrect
    possible_outcome{t,2} = [sound_condition_array(file).VR_sounds{t,5}]; % sound condition
    count_t = count_t+1;
    trial_id = [trial_id,count_t];

end
file_digidata_trial_info(file).estimated_digidata = possible_outcome;
file_digidata_trial_info(file).trial_id = trial_id;


% ex_data = abfload(strcat(digidata_its(file).directory));
% figure(55);clf;
% hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0),'-r'); plot(rescale(ex_data(:,5),-1,0),'-b');plot(rescale(ex_data(:,8),-1,0),'-m');hold off; movegui(gcf,'center');


% count = 0;
% previous_trial = []; % keep track of which gaps were used
% for g = 2:length(difference_between_gaps)
%     
%     if difference_between_gaps(g) > 70*good_dataset_conversion && difference_between_gaps(g) < 101*good_dataset_conversion && ...
%             digidata_its(file).it_gaps_good(big_gaps(g)) > .95*digidata_its(file).sync_sampling_rate && ~ismember(digidata_its(file).locs(big_gaps(g)),previous_trial)
%         count = count+1;
%         possible_outcomes(count).correct = 1;
%         [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g-1)+1) - alignment_info(file).frames_times));
%         possible_outcomes(count).start_trial_digidata_time = digidata_its(file).locs(big_gaps(g-1)+1);
%         if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
%             possible_outcomes(count).frame_start_trial = closest_frame;
%         else
%             possible_outcomes(count).frame_start_trial = nan;
%         end
%         possible_outcomes(count).end_trials_digidata_time = digidata_its(file).locs(big_gaps(g)); %last iteration of current trial during maze (after this is the gap between end of maze and ITI)
%        
%         possible_outcomes(count).start_ITI_digidata_time = digidata_its(file).locs(big_gaps(g)+1); 
%         
%         [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g+1)) - alignment_info(file).frames_times));
%         if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
%             possible_outcomes(count).frame_end_ITI = closest_frame;
%         else
%             possible_outcomes(count).frame_end_ITI = nan;
%         end
%         possible_outcomes(count).end_ITI_digidata_time = digidata_its(file).locs(big_gaps(g+1));
%         previous_trial = [previous_trial, [possible_outcomes(count).start_trial_digidata_time,possible_outcomes(count).end_trials_digidata_time,possible_outcomes(count).start_ITI_digidata_time,possible_outcomes(count).end_ITI_digidata_time]];
%     elseif difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion && digidata_its(file).it_gaps_good(big_gaps(g)) > .8*digidata_its(file).sync_sampling_rate && ...
%             digidata_its(file).it_gaps_good(big_gaps(g)) < .95*digidata_its(file).sync_sampling_rate && ~ismember(digidata_its(file).locs(big_gaps(g)),previous_trial)
%         count = count+1;
%         possible_outcomes(count).correct = 0;
%         possible_outcomes(count).start_trial_digidata_time = digidata_its(file).locs(big_gaps(g-1)+1);
%         [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g-1)+1) - alignment_info(file).frames_times));
%         if val < .02*digidata_its(file).sync_sampling_rate
%             possible_outcomes(count).frame_start_trial = closest_frame;
%         else
%             possible_outcomes(count).frame_start_trial = nan;
%         end
%         possible_outcomes(count).end_trials_digidata_time = digidata_its(file).locs(big_gaps(g));
%         possible_outcomes(count).start_ITI_digidata_time = digidata_its(file).locs(big_gaps(g)+1); 
% 
%         possible_outcomes(count).end_ITI_digidata_time = digidata_its(file).locs(big_gaps(g+1));
%         [val,closest_frame] = min(abs(digidata_its(file).locs(big_gaps(g+1)) - alignment_info(file).frames_times));
%         if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
%             possible_outcomes(count).frame_end_ITI = closest_frame;
%         else
%             possible_outcomes(count).frame_end_ITI = nan;
%         end
%         previous_trial = [previous_trial, [possible_outcomes(count).start_trial_digidata_time,possible_outcomes(count).end_trials_digidata_time,possible_outcomes(count).start_ITI_digidata_time,possible_outcomes(count).end_ITI_digidata_time]];
% 
%         %TO DEAL WITH WEIRD TRIALS!!!!!! where sounds are played in the
%         %wrong place or there are iterations missing from the end of the
%         %trial
%     elseif min(abs(digidata_its(file).locs(big_gaps(g-1)+1) -digidata_time.start_trials_digidata_time)) == 0 && digidata_its(file).locs(big_gaps(g-1)) < digidata_its(file).locs(big_gaps(g-1)+1) ...
%             && digidata_its(file).locs(big_gaps(g-1)+1) < digidata_its(file).locs(big_gaps(g)) && ~ismember(digidata_its(file).locs(big_gaps(g)),previous_trial)
%             [val,trial_id] = min(abs(digidata_time.start_trials_digidata_time - digidata_its(file).locs(big_gaps(g-1)+1)));
%             count = count+1;
%             possible_outcomes(count).weird_trial = count;
%             %assuming I am not loosing as many iterations during the ITI
%             if difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion
%                 possible_outcomes(count).correct = 0;
%             end
% 
%             if difference_between_gaps(g) > 70*good_dataset_conversion && difference_between_gaps(g) < 101*good_dataset_conversion
%                 possible_outcomes(count).correct = 1;
%             end
%             possible_outcomes(count).start_trial_digidata_time = digidata_time.start_trials_digidata_time(trial_id); 
%             [val,closest_frame] = min(abs(digidata_time.start_trials_digidata_time(trial_id) - alignment_info(file).frames_times));
%                 if val < .02*digidata_its(file).sync_sampling_rate
%                     possible_outcomes(count).frame_start_trial = closest_frame;
%                 else
%                     possible_outcomes(count).frame_start_trial = nan;
%                 end
%                 [val,trial_id] =find(digidata_time.end_trials_digidata_time - digidata_time.start_trials_digidata_time(trial_id) >0,1)%min(abs(digidata_time.end_trials_digidata_time - digidata_its(file).locs(big_gaps(g))));
%                 possible_outcomes(count).end_trials_digidata_time = digidata_time.end_trials_digidata_time(trial_id);
%                 [val,trial_id] =find(digidata_time.start_iti_digidata_time - digidata_time.start_trials_digidata_time(trial_id) >0,1);%[val,trial_id] = min(abs(digidata_time.start_iti_digidata_time - digidata_its(file).locs(big_gaps(g)+1)));
%                 possible_outcomes(count).start_ITI_digidata_time = digidata_time.start_iti_digidata_time(trial_id);
% 
%                 [val,trial_id] = min(abs(end_iti_digidata_time - digidata_its(file).locs(big_gaps(g+1))));
%                 possible_outcomes(count).end_ITI_digidata_time = digidata_time.end_iti_digidata_time(trial_id);
%                 [val,closest_frame] = min(abs(digidata_time.end_iti_digidata_time(trial_id)- alignment_info(file).frames_times));
%                 if val < .02*digidata_its(file).sync_sampling_rate %FIND IF THERE IS FRAME 100 or .1ms from this one
%                     possible_outcomes(count).frame_end_ITI = closest_frame;
%                 else
%                     possible_outcomes(count).frame_end_ITI = nan;
%                 end
%                         previous_trial = [previous_trial, [possible_outcomes(count).start_trial_digidata_time,possible_outcomes(count).end_trials_digidata_time,possible_outcomes(count).start_ITI_digidata_time,possible_outcomes(count).end_ITI_digidata_time]];
% 
%     end
%     
% end
%% estimate trials that most closely resemble each other (works for 2 locations)
% 
% % Example data (replace with your actual data)
% trueMatrix = [trial_info.correct; trial_info.condition]'; % Replace with your true trial information matrix
% estimatedMatrix = estimated_trial_info; % Replace with your estimated trial information matrix
% 
% sectionSize = size(estimatedMatrix, 1);
% numSections = size(trueMatrix, 1) - sectionSize + 1;
% mseValues = zeros(numSections, 1);
% 
% for i = 1:numSections
%     trueSection = trueMatrix(i:i+sectionSize-1, :);
%     mseValues(i) = nanmean((trueSection - estimatedMatrix).^2, 'all');
% end
% 
% [minMSE, bestTrueSectionIndex] = min(mseValues);
% 
% file_matching_trials(file,:) = [bestTrueSectionIndex, bestTrueSectionIndex + sectionSize - 1];
% fprintf('Best matching true section of trials: %d to %d\n', bestTrueSectionIndex, bestTrueSectionIndex + sectionSize - 1);
% fprintf('Minimum MSE value: %.4f\n', minMSE);

%%
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
bestTrialIndices = findBestMatchingTrials(trueConditions,estimatedConditions)
file_matching_trials(file,:) = [bestTrialIndices(1), bestTrialIndices(2)];

end
% verify that trials make sense... they go in chronological order
for file = 1:length(digidata_its) 
    if file > 1
        if file_matching_trials(file,1)>=file_matching_trials(file-1,2) == 1
            fprintf('Trial order makes sense!\n'); %\n next line
        else
            fprintf('Error! Trials are out of order! Need to go back and check\n')
        end
    end
end

