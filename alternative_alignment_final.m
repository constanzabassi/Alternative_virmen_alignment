% load data
 mousename = 'HA10-1L';%;
mouse = mousename;
date = '2023-03-31'; %;
server = 'V:';
data_base = 'CBHA10-1L_230331';%;
sync_base_path = [ server '/Connie/RawData/' mousename '/wavesurfer/' date '/'];
virmen_base = [server '/Connie/RawData/' mousename '/virmen/' data_base ];
imaging_base_path=[server '\Connie\RawData\' mousename '\' date '\'];
is_stim_dataset =1; 

data = load(strcat(virmen_base, '.mat'));
dataCell = load(strcat(virmen_base, '_Cell.mat'));

good_dataset = 0;

addpath(genpath('C:\Code\Github\Alternative_virmen_alignment'));

% data = load(strcat(virmen_base, '_1.mat'));
% dataCell = load(strcat(virmen_base, '_Cell_1.mat'));

%% get frame times of all files in this folder
addpath(genpath('C:\Code\Align_signals_imaging'));
galvo_channel = 7;
[alignment_info] = get_frame_times(imaging_base_path, sync_base_path, [], galvo_channel,1,[],[]); %7 is res galvo channel in investigator

%% find sound info for each file
spkr_channel_number = [4,8];
sync_sampling_rate = alignment_info(1).sync_sampling_rate;
distance_between_sounds = 3*sync_sampling_rate ;%min distance between sounds in digidata units in task is ~4 seconds between reward sound and start of trial- (passive at 10k was about 45000)
threshold_spk = 1.4e-4;%0.0001 is what works most of the time!!! increase slightly if not working!!
distance_within_sounds = 0.22*sync_sampling_rate; %for task should be 200
[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task(server,mousename,date,alignment_info,spkr_channel_number,'VR',distance_between_sounds,threshold_spk,distance_within_sounds); %looks at files with VR in the name

%% get trial info using the virmen files!
for tr = 1:length(dataCell.dataCell)
    trial_info(tr).correct = dataCell.dataCell{1,tr}.result.correct;
    trial_info(tr).condition = dataCell.dataCell{1,tr}.maze.condition;
    if is_stim_dataset == 1
        trial_info(tr).is_stim = dataCell.dataCell{1,tr}.maze.is_stim_trial;
    end
end

% use virmen data structure to estimate which iterations match specific
%events (e.g. start_trial_iteration, end_trial_iterations etc.)
%also get the iterations where the mouse crossed 50!
%trial_its has the iteration IDs
[trial_its,trial_its_time] = virmen_it_rough_estimation(data);

%% get digidata iteration locations and difference between them
string = 'VR';
virmen_channel = 6;
digidata_its = get_digidata_iterations(sync_base_path,string, virmen_channel);

%% find its in the data that best match the its for each trial dividing files into trials that match them
[file_estimated_trial_info,file_matching_trials,trial_times,file_digidata_trial_info] = match_trialsperfile(digidata_its, good_dataset, trial_info,sound_condition_array,alignment_info);
% find the start and end trials that are within the imaging frames!// also
% puts trials into context of all other trials (file_digidata_trial_info)
[file_trial_ids,file_digidata_trial_info] = get_trial_ids(file_matching_trials,file_digidata_trial_info,file_estimated_trial_info);
%% start by looking at first trial then shift the real iterations around until they make sense
[possible_alignment] = determine_shift(file_trial_ids,sound_condition_array, trial_its, file_digidata_trial_info,digidata_its,data);
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




