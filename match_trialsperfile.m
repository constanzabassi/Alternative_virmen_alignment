function [file_estimated_trial_info,file_matching_trials,trial_times,file_digidata_trial_info] = match_trialsperfile(digidata_its, good_dataset, trial_info,sound_condition_array,alignment_info)
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
difference_VR = abs(digidata_its(file).locs(big_gaps(g)) - sound_condition_array(file).VR_sounds(:,3));
difference_VR_noabs = (digidata_its(file).locs(big_gaps(g)) - sound_condition_array(file).VR_sounds(:,3));
[~,b]=min((difference_VR)); %finds closest gap
if  any(difference_VR(b) < .2*digidata_its(file).sync_sampling_rate) %might need to change to 0.2 since that might be the greatest difference
%     if difference_VR_noabs(b) < 0 %use sound to say this is where the last gap happens
%         probable_end_trial_gaps(b,1) = sound_condition_array(file).VR_sounds(b,3);
%     else
%         probable_end_trial_gaps(b,1) = [digidata_its(file).locs(big_gaps(g))];
%     end

    probable_end_trial_gaps(b,1) = [digidata_its(file).locs(big_gaps(g))];
    probable_end_trial_gaps(b,2) = [sound_condition_array(file).VR_sounds(b,1)];
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


%%
%start_trials_digidata_time =  digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate)+1);
end_trials_digidata_time = setdiff(digidata_its(file).locs(find(digidata_its(file).it_gaps>.8*digidata_its(file).sync_sampling_rate)),gaps_to_eliminate);
end_iti_digidata_time = setdiff(digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate)),gaps_to_eliminate);%digidata_its(file).locs(find(digidata_its(file).it_gaps >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps < .55*digidata_its(file).sync_sampling_rate));

unpaired = {};
unpaired.iti = setdiff(end_iti_digidata_time,probable_end_iti_gaps(:,1));
unpaired.trial = setdiff(end_trials_digidata_time,probable_end_trial_gaps(:,1));
% pair up the probable end trial and end iti gaps
% pair up so end trial goes first and then end iti
trial_pairs = zeros(max(length(probable_end_trial_gaps),length(probable_end_iti_gaps)),2);
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

start_iti_digidata_time = digidata_its(file).locs(indx+1);
start_trials_digidata_time =  digidata_its(file).locs(indxiti+1);

% start_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps_good>.8*digidata_its(file).sync_sampling_rate)+1);
% start_trials_digidata_time =  digidata_its(file).locs(find(digidata_its(file).it_gaps_good >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps_good < .55*digidata_its(file).sync_sampling_rate)+1);

% build outcome array based on sounds and ITI sounds!
possible_outcome = [];
for t = 1:length(sound_condition_array(file).VR_sounds)
    possible_outcome(t,:) = [sound_condition_array(file).ITI_sounds(t,1),sound_condition_array(file).VR_sounds(t,1)];
end
ex_data = abfload(strcat(digidata_its(file).directory));
figure(55);clf;
hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0),'-r'); plot(rescale(ex_data(:,5),-1,0),'-b');plot(rescale(ex_data(:,8),-1,0),'-m');hold off; movegui(gcf,'center');

ex_data_good = abfload('U:\Connie\RawData\HA10-1L\wavesurfer\2023-03-24\04_VR_8loc_0000.abf');
figure(56);clf;
hold on;plot(ex_data_good(:,7));plot(ex_data_good(:,6)); plot(rescale(ex_data_good(:,4),-1,0),'-r'); plot(rescale(ex_data_good(:,8),-1,0),'-m');hold off; movegui(gcf,'center');

% % Initialize an array to store gaps to eliminate
% gaps_to_eliminate = [];
% % Iterate through gap differences to find gaps too close together
% bad_vals = [];
% for i = 1:length(gap_diffs)-1
%     if gap_diffs(i) < min_threshold
%         % Check if the next gap is far from others
%         if i == 1 || (i > 1 && gap_diffs(i-1) >= min_threshold)
%             gaps_to_eliminate = [gaps_to_eliminate, i+1];
%         else
%             gaps_to_eliminate = [gaps_to_eliminate, i];
%         end
%     end
%     if i >1 && i<length(gap_diffs)-1 && ~ismember(i,gaps_to_eliminate)
%         % Check if the gap differences in front and after the current one are between 3800 and 6000
%         if [(i > 1 && gap_diffs(i-1) >= min_threshold && gap_diffs(i-1) <= max_threshold) && ...
%            (i < length(gap_diffs)-1 && gap_diffs(i+1) >= min_threshold && gap_diffs(i+1) <= max_threshold)]==0
%             bad_vals = [bad_vals,i];
%             %gaps_to_eliminate = [gaps_to_eliminate, i];
%         end
%     end
% end
% dif_bad_vals = diff(bad_vals);
% gaps_to_eliminate = [gaps_to_eliminate, bad_vals(find(dif_bad_vals==1))+2];%plus 2 bc diff twice
% 
% % Iterate through gap differences to find gaps too close together
% % for i = 1:length(gap_diffs)-1
% %     if gap_diffs(i) < min_threshold
% %         % Check if the next gap is far from others
% %         if i == 1 || (i > 1 && gap_diffs(i-1) >= min_threshold)
% %             gaps_to_eliminate = [gaps_to_eliminate, i+1];
% %         else
% %             gaps_to_eliminate = [gaps_to_eliminate, i];
% %         end
% %     end
% % end
% % Eliminate duplicate values from the gaps_to_eliminate array
% gaps_to_eliminate = unique(gaps_to_eliminate);
% big_gaps_all = big_gaps;
% digidata_its(file).locs_corrected = digidata_its(file).locs;
% digidata_its(file).locs_corrected(big_gaps(gaps_to_eliminate)) = [];
% digidata_its(file).it_gaps_good = diff(digidata_its(file).locs_corrected);
% big_gaps(gaps_to_eliminate)=[];
% 
% %try to define trial events using gaps
% end_trials_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps_good>.8*digidata_its(file).sync_sampling_rate));
% start_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps_good>.8*digidata_its(file).sync_sampling_rate)+1);
% % incorrect_trials_digidata_time = find(digidata_its(file).it_gaps_good >.8*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps_good < .95*digidata_its(file).sync_sampling_rate);
% % correct_trials_digidata_time = find(digidata_its(file).it_gaps_good > .95*digidata_its(file).sync_sampling_rate);
% start_trials_digidata_time =  digidata_its(file).locs(find(digidata_its(file).it_gaps_good >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps_good < .55*digidata_its(file).sync_sampling_rate)+1);
% end_iti_digidata_time = digidata_its(file).locs(find(digidata_its(file).it_gaps_good >.25*digidata_its(file).sync_sampling_rate & digidata_its(file).it_gaps_good < .55*digidata_its(file).sync_sampling_rate));
% digidata_time.start_trials_digidata_time = start_trials_digidata_time;
% digidata_time.end_trials_digidata_time = end_trials_digidata_time;
% digidata_time.start_iti_digidata_time = start_iti_digidata_time; 
% digidata_time.end_iti_digidata_time = end_iti_digidata_time;
% trial_times(file) = digidata_time;
% 
% %write code to say if this big gap is followed by another big gap at this time point then!
% difference_between_gaps = diff(big_gaps);
% 
% possible_outcomes = {};
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
% %check to make sure we are not skipping any trials!!
% ex_start = start_trials_digidata_time;
% ex_start(end) = []; %delete bc it would not be in the difference matrix
% no_skips = ismember(ex_start,[possible_outcomes.start_trial_digidata_time]);
% if sum(find(no_skips == 0))>1
%     fprintf('Something is wrong, might be skipping trials\n');
% else
%     fprintf('Trials seem fine\n');
% end
% 
% ex_data = abfload(strcat(digidata_its(file).directory));
% if ~isempty(gaps_to_eliminate)
%     
% figure(file);clf;title(strcat('File ',num2str(file),' has weird gaps (deleted in red)'))
% hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(digidata_its(file).locs(big_gaps_all),0,'*c');plot(digidata_its(file).locs(big_gaps_all(gaps_to_eliminate)),0,'*r');hold off; movegui(gcf,'center');
% else
%     figure(file);clf;
%     title(strcat(num2str(file),' file has no weird gaps'))
%     hold on;plot(ex_data(:,7));plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(digidata_its(file).locs(big_gaps_all),0,'*c');hold off; movegui(gcf,'center');
% end
% % combine sound condition and trial outcome based on digidata time!
% % finds difference between start trial and sound offset
% estimated_trial_info = [];weird_trials = []; count_t=0;trial_id = [];
% for t = 1:length(possible_outcomes)
%     %[val,closest_sound] = min(abs(possible_outcomes(t).digidata_time - sound_condition_array(file).file(:,3)));
%     if ~isnan(possible_outcomes(t).frame_start_trial)== 1%if it's 1 then the trial might be too short to look at
%         count_t = count_t+1;
%         [closest_sound] = find(possible_outcomes(t).end_trials_digidata_time - sound_condition_array(file).file(:,3)<0,1)-1; %find first value where it becomes negative and use the one before it
%         estimated_trial_info(count_t,:) = [possible_outcomes(t).correct;sound_condition_array(file).file(closest_sound,1)];
%         trial_id(count_t) = t;
%         if sound_condition_array(file).file(closest_sound,2)< possible_outcomes(t).start_trial_digidata_time  %if the sound plays during the ITI make it a NaN
%                 weird_trials = [weird_trials,count_t];
%                 estimated_trial_info(count_t,:) = [possible_outcomes(t).correct;nan];
%         end
%     end
% end
% file_estimated_trial_info(file).estimated_trial_info = estimated_trial_info;
% file_estimated_trial_info(file).trial_id = trial_id;
% file_estimated_trial_info(file).weird_trials = weird_trials;
% file_digidata_trial_info(file).estimated_digidata = possible_outcomes;
% %figure(); hold on;plot(ex_data(:,6));plot(possible_correct_trial,ex_data(possible_correct_trial,6),'*c');plot(possible_incorrect_trial,ex_data(possible_incorrect_trial,6),'*r');hold off
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
