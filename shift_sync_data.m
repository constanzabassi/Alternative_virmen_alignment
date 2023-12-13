function [virmen_it,trial_its] = shift_sync_data(data,file_trial_ids,digidata_its,file_estimated_trial_info,sound_condition_array,task_info)
[trial_its,trial_its_time] = virmen_it_rough_estimation(data); 

for file = 1:length(digidata_its)
%determine the timing of each virmen iteration in the session
iterations_in_time = data.data(1,:).*(86400);
iterations_in_time = iterations_in_time.*digidata_its(file).sync_sampling_rate;
it_ids = 1:length(data.data(1,:));
mean_freq = round(mode(diff(iterations_in_time))); 

start_trial_number = file_trial_ids(file,1);%+file_trial_ids(file,3)-2; 
end_trial_number = file_trial_ids(file,2);

%initialize variables
possible_it_times = [];
possible_iterations = [];
possible_it_locs =[];
last_iteration_in_file = [];
shift = [];pos_peak_id =[];pos_peak_pos=[];
% case 1- I have a positive iteration! use this to determine shift and to
% number iterations accordingly
if ~isempty(digidata_its(file).pos_loc)
    pos_peak_id = digidata_its(file).pos_pks(1)*10e4;
    pos_peak_id = round(pos_peak_id);
    
    if pos_peak_id(1) == 1000000 %very first file has large positive peak (not actually an iteration of virmen have to add 1!)
        % make sure next iteration is close/ makes sense
        if digidata_its(file).locs(find(digidata_its(file).locs == digidata_its(file).pos_loc(1))+1) - digidata_its(file).locs(find(digidata_its(file).locs == digidata_its(file).pos_loc(1)))< mean_freq*2
            pos_peak_pos = digidata_its(file).locs(find(digidata_its(file).locs == digidata_its(file).pos_loc(1))+1);
        else
            pos_peak_pos = digidata_its(file).pos_loc(1) + mean_freq;
        end
        
        shift = pos_peak_pos;% need to test vs pos_peak_pos  pos_peak_pos;%
        possible_it_times = iterations_in_time+shift-iterations_in_time(1); %assign iterations times with added shift
        %iteration ids of iterations withing imaging limits
        possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits
        
        %test distance of unfinished sounds and last iteration of trial to
        %see if they are close together otherwise add a shift of 1
        %iteration
        test_iterations = trial_its.end_trial_its(start_trial_number:end_trial_number);
        all_sound_trials = 1+file_trial_ids(file,3):length(test_iterations)+file_trial_ids(file,3);
        excluded_trials = find([sound_condition_array(file).VR_sounds{:,6}]~=1); %only include unifinished sounds in analysis
        [sound_trials,ind] = setdiff(all_sound_trials,excluded_trials);
        difference_it_sound = [([sound_condition_array(file).VR_sounds{sound_trials,3}]-0.005*digidata_its(file).sync_sampling_rate) - possible_it_locs(test_iterations(ind)-possible_iterations(1)+1)];
        small_shift = round(mean(difference_it_sound(find(difference_it_sound < mean_freq*3 & difference_it_sound > 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift
        small_shift_neg = round(mean(difference_it_sound(find(difference_it_sound > -mean_freq*5 & difference_it_sound < 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift

        if ~isnan(small_shift) && isnan(small_shift_neg)
            new_shift = shift+small_shift;
        elseif ~isnan(small_shift_neg) && isnan(small_shift)
            new_shift = shift+small_shift_neg;
        else
            new_shift = shift;%
        end


        possible_it_times = iterations_in_time+new_shift-iterations_in_time(1); %assign iterations times with added shift
        %iteration ids of iterations withing imaging limits
        possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits
        
    else
        pos_peak_pos = digidata_its(file).pos_loc(1);%digidata_its(file).locs(find(digidata_its(file).locs == digidata_its(file).pos_loc(1)))+1);%digidata_its(file).pos_loc(1);
        shift = pos_peak_pos; %peak position
        possible_it_times = iterations_in_time+shift-iterations_in_time(pos_peak_id-1); %assign iterations times with added shift
        %iteration ids of iterations withing imaging limits
        possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits
        
        %test distance of unfinished sounds and last iteration of trial to
        %see if they are close together otherwise add a shift of 1
        %iteration
        test_iterations = trial_its.end_trial_its(start_trial_number:end_trial_number);
        all_sound_trials = 1+file_trial_ids(file,3):length(test_iterations)+file_trial_ids(file,3);
        excluded_trials = find([sound_condition_array(file).VR_sounds{:,6}]~=1); %only include unifinished sounds in analysis
        [sound_trials,ind] = setdiff(all_sound_trials,excluded_trials);
        difference_it_sound = [([sound_condition_array(file).VR_sounds{sound_trials,3}]-0.005*digidata_its(file).sync_sampling_rate) - possible_it_locs(test_iterations(ind)-possible_iterations(1)+1)];
        % sound always ends ~3ms after last iteration in the trial
        small_shift = round(mean(difference_it_sound(find(difference_it_sound < mean_freq*3 & difference_it_sound > 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift
        small_shift_neg = round(mean(difference_it_sound(find(difference_it_sound > -mean_freq*5 & difference_it_sound < 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift

        if ~isnan(small_shift) && isnan(small_shift_neg)
            new_shift = shift+small_shift;
        elseif ~isnan(small_shift_neg) && isnan(small_shift)
            new_shift = shift+small_shift_neg;
        else
            new_shift = shift;%
        end


        possible_it_times = iterations_in_time+new_shift-iterations_in_time(pos_peak_id-1); %assign iterations times with added shift
        %iteration ids of iterations withing imaging limits
        possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
        possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits
    end
else
    % first say that the first iteration time in the first full trial is in the
    % position of the first iteration in the trua dataset (virmen data)
    shift = file_estimated_trial_info(file).start_trials_digidata_time(file_trial_ids(file,3));%(start_trial_number); 
    
    possible_it_times = iterations_in_time+shift-iterations_in_time(trial_its.start_trial_its(start_trial_number));
    %start by assuming that the estimated digidata time alings perfectly with
    %true data- here is where I figure out how to shift over time
    possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
    possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits
    
    %test distance of unfinished sounds and last iteration of trial to
    %see if they are close together otherwise add a shift of 1
    %iteration
    test_iterations = trial_its.end_trial_its(start_trial_number:end_trial_number);
    all_sound_trials = 1+file_trial_ids(file,3):length(test_iterations)+file_trial_ids(file,3);
    excluded_trials = find([sound_condition_array(file).VR_sounds{:,6}]~=1); %only include unifinished sounds in analysis
    [sound_trials,ind] = setdiff(all_sound_trials,excluded_trials);
    difference_it_sound = [([sound_condition_array(file).VR_sounds{sound_trials,3}]-0.005*digidata_its(file).sync_sampling_rate) - possible_it_locs(test_iterations(ind)-possible_iterations(1)+1)];
    small_shift = round(mean(difference_it_sound(find(difference_it_sound < mean_freq*3 & difference_it_sound > 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift
    small_shift_neg = round(mean(difference_it_sound(find(difference_it_sound > -mean_freq*5 & difference_it_sound < 0)))/(mean_freq/2))*(mean_freq/2); %find closest ones and determine if there needs to be another small shift

    if ~isnan(small_shift) && isnan(small_shift_neg)
        new_shift = shift+small_shift;
    elseif ~isnan(small_shift_neg) && isnan(small_shift)
        new_shift = shift+small_shift_neg;
    else
        new_shift = shift;
    end

    possible_it_times = iterations_in_time+new_shift-iterations_in_time(trial_its.start_trial_its(start_trial_number));
    possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames// this is iterations ids
    possible_it_locs = possible_it_times(possible_iterations); %locations of iterations within limits

end

%final test- see if there are about 7 iterations from tiny gap before sound
%onset
sound_onsets_speakers = [sound_condition_array(file).VR_sounds{:,2}]; 
nan_ind = find(isnan(sound_onsets_speakers));
sound_onsets_speakers = sound_onsets_speakers(~isnan(sound_onsets_speakers));
sound_onsets_iterations = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > trial_its.start_trial_its(start_trial_number) & trial_its.sound_trigger_its <trial_its.end_iti_its(end_trial_number)))+6; %sound happens within 7 iterations
possible_sound_onsets = possible_it_times(sound_onsets_iterations);

all_differences = [];
for s = 1:length(possible_sound_onsets)
    difference_sounds = min(abs(possible_sound_onsets(s) - sound_onsets_speakers));
    all_differences = [all_differences,difference_sounds];
%     if difference_sounds < 0.012 * digidata_its(file).sync_sampling_rate %if they are within 12ms of each other
%         fprintf('Sound distances make sense!\n');
%     else
%         fprintf('Sound distances do not make sense!\n');
%     end
end

%save average distance from trigger iteration right after the gap (trigger
%+1)
sound_onsets_iterations2 = trial_its.sound_trigger_its(find(trial_its.sound_trigger_its > trial_its.start_trial_its(start_trial_number) & trial_its.sound_trigger_its <trial_its.end_iti_its(end_trial_number)))+0; %sound happens within 7 iterations
possible_sound_onsets2 = possible_it_times(sound_onsets_iterations2);
% mean_sound_distance = mean(sound_onsets_speakers(all_sound_trials)-possible_sound_onsets2);% works only if there are no NANs sounds           
mean_sound_distance = [];
for s = 1:length(possible_sound_onsets2)
    difference_sounds_it = min(abs(possible_sound_onsets2(s) - sound_onsets_speakers));
    mean_sound_distance = [mean_sound_distance,difference_sounds_it]; %gives large values for trials without sounds
end
mean_sound_distance = mean_sound_distance(find(mean_sound_distance <.2*digidata_its(file).sync_sampling_rate));
mean_sound_distance = mean(mean_sound_distance);

figure(998);clf; 
hold on
title(strcat('Sound onset verification less than 100ms apart -file # ', num2str(file)));
histogram(all_differences(find(all_differences< 0.1* digidata_its(file).sync_sampling_rate)),'BinWidth',2); %
xline(mean_freq,'-r')
xlabel('Distance between iteration at onset and sound onset in ms')
ylabel('Number of sound onsets')
hold off


ex_data = abfload(strcat(digidata_its(file).directory));
figure(999);clf; 
title(strcat('Shifted data file # ', num2str(file)));
hold on; aa = plot(ex_data(:,task_info.channel_number(1)));bb = plot(ex_data(:,task_info.channel_number(2)),'-k');  cc = plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-b');dd = plot(rescale(ex_data(:,task_info.channel_number(4)),-1,0),'-m');a = plot(possible_it_locs,0,'*c');
legend([aa bb cc dd  a(1)],'Imaging frames','Virmen its','Speaker 1','Speaker 2', 'Estimated iteration times')
if ~isempty(pos_peak_id) && pos_peak_id-possible_iterations(1) > 0 && pos_peak_id-possible_iterations(1) < length(possible_it_locs)
    if pos_peak_id(1) == 1000000
        plot(possible_it_locs(1),0,'*r');
    else
        plot(possible_it_locs(pos_peak_id-possible_iterations(1)),0,'*r');
    end
end
hold off
pause
%hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(possible_it_locs,0,'*c');hold off; movegui(gcf,'center');

virmen_it(file).it_times = possible_it_locs;
virmen_it(file).actual_it_values = possible_iterations; %ids
virmen_it(file).shift = new_shift;
virmen_it(file).difference = all_differences;
virmen_it(file).directory = digidata_its(file).directory;
virmen_it(file).start_trial_number = start_trial_number;
virmen_it(file).end_trial_number = end_trial_number;

virmen_it(file).file_trial_id_start = file_trial_ids(file,3);
virmen_it(file).file_trial_id_end = file_trial_ids(file,4);

virmen_it(file).mean_sound_distance = mean_sound_distance;


end

% figure(999);clf; hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(rescale(ex_data(:,5),-1,0));plot(possible_it_locs(test_iterations-possible_iterations(1)+1),0,'*c');hold off; movegui(gcf,'center');
% %figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');
% figure(999);clf; hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(rescale(ex_data(:,8),-1,0));plot(rescale(ex_data(:,5),-1,0));plot(possible_it_locs(test_iterations-possible_iterations(1)+1),0,'*c');hold off; movegui(gcf,'center');
end