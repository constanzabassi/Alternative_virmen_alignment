% load data
mousename = 'HA11-1R';%;
mouse = mousename;
date = '2023-04-07'; %;
server = 'V:';
data_base = 'CBHA11-1R_230407';%;
sync_base_path = [ server '/Connie/RawData/' mousename '/wavesurfer/' date '/'];
virmen_base = [server '/Connie/RawData/' mousename '/virmen/' data_base ];
imaging_base_path=[server '\Connie\RawData\' mousename '\' date '\'];
is_stim_dataset = 1; 
% give data inputs!
galvo_channel = 7;
good_dataset = 0;
addpath(genpath('C:\Code\Github\Alternative_virmen_alignment'));
cd('C:\Code\Github\Alternative_virmen_alignment');
addpath(genpath('C:\Code\Align_signals_imaging'));

%% load virmen data

data = load(strcat(virmen_base, '.mat'));
dataCell = load(strcat(virmen_base, '_Cell.mat'));


% data = load(strcat(virmen_base, '_1.mat'));
% dataCell = load(strcat(virmen_base, '_Cell_1.mat'));

%% get frame times of all files in this folder
mkdir(strcat(server,'\Connie\ProcessedData\',num2str(mouse),'\',num2str(date)))
cd(strcat(server,'\Connie\ProcessedData\',num2str(mouse),'\',num2str(date)));
if isfile("alignment_info.mat")
    load("alignment_info.mat");
else
    [alignment_info] = get_frame_times(imaging_base_path, sync_base_path, [], galvo_channel,1,[],[]); %7 is res galvo channel in investigator
end
save ('alignment_info','alignment_info');
%% find sound info for each file
sound_info = {};
sound_info.spkr_channel_number = [4,8];%[4,5,8];
sound_info.speaker_ids = [1,2];%[1,2,4]; 
sound_info.mult_spkr = 0; %if multiple speakers are used in a single trial (8 locs)
sound_info.detection_threshold = 0.4;%for 1k (0.45)between 0.4 and 0.5 (0.5 gets rid of more noise) - for some 10k 0.8 (one file #8 in HA10-1L\2023-03-24)

sound_info.sync_sampling_rate = alignment_info(1).sync_sampling_rate;
sound_info.distance_between_sounds = 2*sound_info.sync_sampling_rate ;%min distance between sounds in digidata units in task is ~4 seconds between reward sound and start of trial- (passive at 10k was about 45000)
sound_info.distance_within_sounds = 0.2*sound_info.sync_sampling_rate; %for task should be 200
sound_info.sound_duration = [0.99*sound_info.sync_sampling_rate,1.1*sound_info.sync_sampling_rate];
sound_info.correct = .250; %correct_trial_ITI_length in seconds
sound_info.incorrect = .40; %incorrect_trial_ITI_length in seconds
sound_info.smoothing_factor = 15; %almost always 15 sometimes 20


[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task_new(server,mousename,date,alignment_info,'VR',sound_info);

%[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task_new(server,mousename,date,alignment_info,spkr_channel_number,'VR',detection_threshold,distance_between_sounds,distance_within_sounds,sound_duration,correct_trial_ITI_length,incorrect_trial_ITI_length,mult_spkr,smoothing_factor) 
%pc=1 if windows, any other number if mac

%% get trial info using the virmen files!
for tr = 1:length(dataCell.dataCell)
    trial_info(tr).correct = dataCell.dataCell{1,tr}.result.correct;
    trial_info(tr).condition = dataCell.dataCell{1,tr}.maze.condition;
    if is_stim_dataset == 1
        trial_info(tr).is_stim = dataCell.dataCell{1,tr}.maze.is_stim_trial;
    end
end

%% get digidata iteration locations and difference between them
string = 'VR';
virmen_channel = 6;
digidata_its = get_digidata_iterations(sync_base_path,string, virmen_channel);

%% find its in the data that best match the its for each trial dividing files into trials that match them
task_info.correct = 3; %correct ITI time in sec (distance between end trial and end ITI is ITI +1 sec due to virmen bug)
task_info.incorrect = 5; %incorrect ITI time in sec
task_info.min = 3.5;% minimum time in sec to complete a trial (without ITI)
task_info.channel_number = [galvo_channel,virmen_channel,sound_info.spkr_channel_number];

[file_estimated_trial_info,file_matching_trials] = match_trialsperfile(digidata_its, trial_info,sound_condition_array,task_info);
%[file_estimated_trial_info,file_matching_trials,trial_times,file_digidata_trial_info] = match_trialsperfile(digidata_its, good_dataset, trial_info,sound_condition_array,alignment_info);

%determine closest frames for each trial event
%file_estimated_trial_info = determine_close_frames (alignment_info,file_estimated_trial_info);

% find the start and end trials that are within the imaging frames!// also
% puts trials into context of all other trials (file_digidata_trial_info)
[file_trial_ids,file_digidata_trial_info] = get_trial_ids(file_matching_trials,file_estimated_trial_info,alignment_info,sync_base_path,task_info);
%% shift iterations in time until they match positive peaks or first trial iteration
[virmen_it,trial_its] = shift_sync_data(data,file_trial_ids,digidata_its,file_estimated_trial_info,sound_condition_array,task_info);


%[possible_alignment] = determine_shift(file_trial_ids,sound_condition_array, trial_its, file_digidata_trial_info,digidata_its,data);
%figure(); hold on;plot(ex_data(:,6));plot(ex_data(:,7)); plot(rescale(ex_data(:,4),-1,0));plot(possible_alignment(1).it_times,-.5,'*c');plot(possible_alignmen(1).sound_onsets,-.5,'*r');hold off; movegui(gcf,'center');
%% try to estimate iterations based on positive peaks if they exist
% [positive_peaks,positive_locs] = findpeaks(ex_data(:,virmen_channel),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5);
% positive_peaks = round(positive_peaks,1)*1e5;
% data.data(1,positive_peaks(1))

%% TIME CALCULATIONS
% clear trial_time trial_ex
% 
% for t = 1:size(start_trial_its, 2)
%     t
%     trial_duration = data.data(1, end_iti_its(t)) - data.data(1, 1); % Calculate the duration -multiply by * 86400 gives seconds!   
%     % Format the trial duration depending on whether it's over an hour
%     trial_time_fromstart(t, :) = trial_duration* 86400;
%     if t < size(start_trial_its, 2)
%         trial_time(t,:) = [data.data(1, end_iti_its(t+1)) - data.data(1, start_trial_its(t))] *86400;
%     else
%         trial_time(t,:) = nan;
%     end
% end
% 
% % calculate time differences between consecutive abf files!
% for f = 1:size(sync_dir, 1)
%     session_duration = datestr(sync_dir(f).date, 'HH.MM.SS');    
%     file_time(f, :) = session_duration; %str2double((session_duration));
% 
% end
% 
% % Convert time values to serial date numbers
% serial_numbers = datenum(file_time, 'HH.MM.SS');
% 
% % Calculate time differences in seconds
% file_time_diff_seconds = diff(serial_numbers) * 86400; % Convert days to seconds (1 day = 24*60*60 seconds)




