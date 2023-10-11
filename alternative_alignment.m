% load data
mousename = 'HE4-1L1R';%'HA10-1L';%;
mouse = mousename;
date = '2023-08-24';%'2023-03-30'; %;
server = 'U:';
data_base = 'CBHE4-1L1R_230824';%'CBHA10-1L_230330';%;
sync_base_path = [ server '/Connie/RawData/' mousename '/wavesurfer/' date '/'];
virmen_base = [server '/Connie/RawData/' mousename '/virmen/' data_base ];
imaging_base_path=[server '\Connie\RawData\' mousename '\' date '\'];
is_stim_dataset =0; 

data = load(strcat(virmen_base, '.mat'));
dataCell = load(strcat(virmen_base, '_Cell.mat'));

good_dataset = 1;

% data = load(strcat(virmen_base, '_1.mat'));
% dataCell = load(strcat(virmen_base, '_Cell_1.mat'));

%% get frame times
addpath 'C:\Users\RUNYAN1\Documents\Align_signals_imaging';
[alignment_info] = get_frame_times(imaging_base_path, sync_base_path, [], 7,1,[],[]); %7 is res galvo channel in investigator
%% get trial info
for tr = 1:length(dataCell.dataCell)
    trial_info(tr).correct = dataCell.dataCell{1,tr}.result.correct;
    trial_info(tr).condition = dataCell.dataCell{1,tr}.maze.condition;
    if is_stim_dataset == 1
        trial_info(tr).is_stim = dataCell.dataCell{1,tr}.maze.is_stim_trial;
    end
end

%% look at digidata data and do super rough alignment
virmen_channel = 6;
[ex_data, sync_sampling_interval] = abfload(strcat(sync_base_path,'02_VR_2locswstim_0000.abf')); %'06_VR_2locwstim_0001.abf
sync_sampling_rate = 1/sync_sampling_interval*1e6;
[pks,locs] = findpeaks(abs(ex_data(:,virmen_channel)),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5);
it_gaps = diff(locs);

%look at sound onset first
threshold = 350; %lowest value for this mouse is ~454 in 500 maze- distance between sounds on different trials??
temp_onset = sort([find(round(floor(data.data(3,:))) == 50),find(round(floor(data.data(3,:))) == 51),find(round(floor(data.data(3,:))) == 52)]);
%find ones that are too close together and delete them
temp = find(diff(temp_onset)< threshold);
temp_onset(temp+1) = [];
sound_onset_its = temp_onset;
%figure();hold on;plot(ex_data(:,6));plot(rescale(ex_data(:,4),-1,0));plot(virmen_aq(1).it_times(sound_onset_its(8:28)-it_num_first_one_in_file),-.5,'*c');plot(ex_data(:,7));hold off

%look at iti to determine isITI and also start of trial
start_iti_its = find(diff(data.data(9,:))>0);
start_iti_its = start_iti_its +1; 
end_iti_its = find(diff(data.data(9,:))<0);
start_trial_its = end_iti_its + 1;
start_trial_its = sort([start_trial_its,1]); %making sure I have an it for the very first trial
end_trial_its = start_iti_its - 1;

%get complete trials
min_trials = min([length(end_trial_its),length(start_trial_its)]);
if length(start_trial_its) > min_trials
    start_trial_its(end) = [];
end
%% looking at iteration timing info! try to tie it back to the big sterotyped gaps that I see!
it_times = data.data(1,:).* 86400; %convert to seconds
it_time_gaps = diff(it_times);

%look at iti to determine isITI and also start of trial
end_trial_its_time = it_times(find(it_time_gaps>.8));
start_iti_its_time = it_times(find(it_time_gaps>.8)+1);
incorrect_its_digidata_time = it_times(find(it_time_gaps >.8 & it_time_gaps < .95));
correct_its_its_time = it_times(find(it_time_gaps > .95));
temp_its = find(it_time_gaps > .25 & it_time_gaps < .51)+1;
start_trial_its_time = it_times([1,temp_its]); %adding one for very first trial

%% TIME CALCULATIONS
clear trial_time trial_ex trial_time_fromstart

for t = 1:size(start_trial_its, 2)
    t
    trial_duration = data.data(1, end_iti_its(t)) - data.data(1, 1); % Calculate the duration -multiply by * 86400 gives seconds!   
    % Format the trial duration depending on whether it's over an hour
    trial_time_fromstart(t, :) = trial_duration* 86400;
    if t < size(start_trial_its, 2)
        trial_time(t,:) = [data.data(1, end_iti_its(t+1)) - data.data(1, start_trial_its(t))] *86400;
    else
        trial_time(t,:) = nan;
    end
end

% calculate time differences between consecutive abf files!
for f = 1:size(sync_dir, 1)
    session_duration = datestr(sync_dir(f).date, 'HH.MM.SS');    
    file_time(f, :) = session_duration; %str2double((session_duration));

end

% Convert time values to serial date numbers
serial_numbers = datenum(file_time, 'HH.MM.SS');

% Calculate time differences in seconds
file_time_diff_seconds = diff(serial_numbers) * 86400; % Convert days to seconds (1 day = 24*60*60 seconds)

%% find sound info for each file
spkr_channel_number = [4,8];
distance_between_sounds = 3*sync_sampling_rate ;%min distance between sounds in digidata units in task is ~4 seconds between reward sound and start of trial- (passive at 10k was about 45000)
[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task(server,mousename,date,alignment_info,spkr_channel_number,'VR',distance_between_sounds); %looks at files with VR in the name

%% find its in the data that best match the its for each trial dividing files into trials that match them
file = 2;
end_trial_digidata_time = locs(find(it_gaps>.8*sync_sampling_rate));
start_iti_digidata_time = locs(find(it_gaps>.8*sync_sampling_rate)+1);
incorrect_trials_digidata_time = find(it_gaps >.8*sync_sampling_rate & it_gaps < .95*sync_sampling_rate);
correct_trials_digidata_time = find(it_gaps > .95*sync_sampling_rate);
start_trials_digidata_time =  locs(find(it_gaps >.25*sync_sampling_rate & it_gaps < .55*sync_sampling_rate)+1);
big_gaps = find(it_gaps > .25*sync_sampling_rate); 
% figure(); hold on; plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(locs((big_gaps)),0,'*c'); plot(ex_data(:,7)); hold off
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
possible_correct_trial = [];
possible_incorrect_trial = [];
possible_outcomes = {};
count = 0;

for g = 1:length(difference_between_gaps)
    if difference_between_gaps(g) > 70*good_dataset_conversion && difference_between_gaps(g) < 101*good_dataset_conversion && it_gaps(big_gaps(g)) > .95*sync_sampling_rate
        possible_correct_trial = [possible_correct_trial,locs(big_gaps(g))];
        [val,closest_frame] = min(abs(locs(big_gaps(g)) - alignment_info(file).frames_times));
        count = count+1;
        possible_outcomes(count).correct = 1;
        possible_outcomes(count).digidata_time = locs(big_gaps(g));
        if val < .02*sync_sampling_rate
            possible_outcomes(count).frame = closest_frame;
        else
            possible_outcomes(count).frame = nan;
        end
    elseif difference_between_gaps(g) > 120*good_dataset_conversion && difference_between_gaps(g) < 160*good_dataset_conversion && it_gaps(big_gaps(g)) > .8*sync_sampling_rate && it_gaps(big_gaps(g)) < .95*sync_sampling_rate
        possible_incorrect_trial = [possible_incorrect_trial,locs(big_gaps(g))];
        [val,closest_frame] = min(abs(locs(big_gaps(g)) - alignment_info(file).frames_times));
        count = count+1;
        possible_outcomes(count).correct = 0;
        possible_outcomes(count).digidata_time = locs(big_gaps(g));
        if val < .02*sync_sampling_rate
            possible_outcomes(count).frame = closest_frame;
        else
            possible_outcomes(count).frame = nan;
        end
    end
    
end
% combine sound condition and trial outcome based on digidata time!
estimated_trial_info = [];
for t = 1:length(possible_outcomes)
    [val,closest_sound] = min(abs(possible_outcomes(t).digidata_time - sound_condition_array(file).file(:,3)));
    estimated_trial_info(t,:) = [possible_outcomes(t).correct;sound_condition_array(file).file(closest_sound,1)];
end
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
    mseValues(i) = mean((trueSection - estimatedMatrix).^2, 'all');
end

[minMSE, bestTrueSectionIndex] = min(mseValues);

fprintf('Best matching true section of trials: %d to %d\n', bestTrueSectionIndex, bestTrueSectionIndex + sectionSize - 1);
fprintf('Minimum MSE value: %.4f\n', minMSE);

%% start by looking at first trial then shift the real iterations around until they make sense
start_trial_number = 64; end_trial_number = 86;
%start_trial_its is the iti #
%possible iterations within file
possible_iterations = start_trial_its(start_trial_number): end_iti_its(end_trial_number);

% first say that the first iteration time in the first full trial is in the
% position of the first iteration in the trua dataset (virmen data)
possible_digidata_time_start = start_trials_digidata_time(2);%(start_trial_number); 
%it_time_gaps
%use time difference and assume first time difference to be exactly zero
it_times_this_file = [data.data(1,possible_iterations)-data.data(1,possible_iterations(1))].* 86400;% difference with first iteration for each iteration in seconds!
it_times_this_file = it_times_this_file*sync_sampling_rate;

%start by assuming that the estimated digidata time alings perfectly with
%true data- here is where I figure out how to shift over time
possible_it_times = it_times_this_file + possible_digidata_time_start;
figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');



%% try to estimate iterations based on positive peaks if they exist
% [positive_peaks,positive_locs] = findpeaks(ex_data(:,virmen_channel),'MINPEAKHEIGHT',0.09,'MinPeakDistance',5);
% positive_peaks = round(positive_peaks,1)*1e5;
% data.data(1,positive_peaks(1))





