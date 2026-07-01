%% provide all inputs
info.mousename = 'KN8-3L';%;
info.mouse = info.mousename;
info.date = '2026-06-24'; %;
info.server = 'W:'; %/Volumes/Runyan5
info.mouse_date = 'KN8-3L/2026-06-24';
runyan5 = "V:";
runyan4 = 'W:';
data_base = 'AGKN-8-3L_260624';%;
info.sync_base_path = [ info.server '/Connie/RawData/' info.mousename '/wavesurfer/' info.date '/'];
% info.virmen_base = [info.server '/Connie/RawData/' info.mousename '/virmen/' data_base ];
info.imaging_base_path=[info.server '/Connie/RawData/' info.mousename '/' info.date '/'];
info.save_path = [info.server '/Connie/ProcessedData/' info.mousename '/' info.date '/passive/'];
info.processed_path = [info.server '/Connie/ProcessedData/' info.mousename '/' info.date '/'];
info.is_stim_dataset = 1; 
% give data inputs!
info.galvo_channel = 6;
info.virmen_channel = 5;

sound_info = {}; 
sound_info.spkr_channel_number = [4,7,8];%[4,5,8];
sound_info.speaker_ids = [1,3,2];%[1,2,4]; 
sound_info.mult_spkr = 1; %if multiple speakers are used in a single trial (8 locs)
%load conditions per speaker in runyan 5 info.server
load(strcat(runyan5,'/Connie/condition_per_speaker'));
%load('/Volumes/Runyan5/Connie/condition_per_speaker.mat');
sound_info.condition_per_speaker = conditions_per_speaker;
sound_info.info = info;


code_folder = uigetdir; 
addpath(genpath(code_folder));


cd(strcat(num2str(ss),'\Connie\ProcessedData\',num2str(info.mouse_date),'\VR\'));
load('alignment_variables.mat','info','sound_info');
info.server = ss;

info.sync_base_path = [ info.server '/Connie/RawData/' info.mousename '/wavesurfer/' info.date '/'];
% info.virmen_base = [info.server '/Connie/RawData/' info.mousename '/virmen/' data_base ];
info.imaging_base_path=[info.server '/Connie/RawData/' info.mousename '/' info.date '/'];
info.save_path = [info.server '/Connie/ProcessedData/' info.mousename '/' info.date '/passive/'];
info.processed_path = [info.server '/Connie/ProcessedData/' info.mousename '/' info.date '/'];

%% LOAD ALIGNMENT DATA!
% mkdir(strcat(info.server,'/Connie/ProcessedData/',num2str(info.mouse),'/',num2str(info.date)))
cd(strcat(info.server,'/Connie/ProcessedData/',num2str(info.mouse),'/',num2str(info.date)));
if isfile("alignment_info.mat")
    load("alignment_info.mat");
end
%% FIND SOUNDS!

%passive data sounds
info.save_path = [info.server '/Connie/ProcessedData/' info.mousename '/' info.date '/passive/'];

sound_info.sync_sampling_rate = alignment_info(1).sync_sampling_rate;
sound_info.distance_between_sounds = 2*sound_info.sync_sampling_rate ;%min distance between sounds in digidata units in task is ~4 seconds between reward sound and start of trial- (passive at 10k was about 45000)
d1 = datetime(date);
d2 = datetime('2023-07-03'); %when sensors where changed in the investigator
if d1 < d2
    sound_info.distance_within_sounds = 0.05*sound_info.sync_sampling_rate; %for task should be 200/sometimes 50
else
    sound_info.distance_within_sounds = 0.2*sound_info.sync_sampling_rate; %for task should be 200/sometimes 50
end
sound_info.sound_duration = 1*sound_info.sync_sampling_rate;%[0.99*sound_info.sync_sampling_rate,1.1*sound_info.sync_sampling_rate];
sound_info.correct = []; 
sound_info.incorrect = []; 
sound_info.smoothing_factor = 15; %almost always 15 sometimes 20

sound_info.unique_detection_threshold = [];%list specific file and threshold wanted [file#1,threshold1; file#2,threshold2]
sound_info.detection_threshold = 1.5;%for 1k (0.45)between 0.4 and 0.5 (0.5 gets rid of more noise) - for some 10k 0.8 (one file #8 in HA10-1L\2023-03-24)

[sound_st, sound_trials, sound_condition_array] = find_spkr_output_task_simple(info,alignment_info,'passive',sound_info);

%% FIND FRAMES FOR EACH SOUND (AND BINARIZE)
[passive_frames,new_sound_st,sounds_per_file] = binarize_passive_sounds(sound_st,sound_info,alignment_info);
%passive_frames.corr_frames = [onset offset] of each sound in terms of all other frames in imaging session
%passive_frames.frames = [onset offset] of each sound in terms of frames within file
%passive_frames.trial_num puts repeats together into same trial
%passive_frames.og_trial_num puts repeats together into same trial in terms of frames within file
%passive_frames.condition is the true sound condition

%% SAVE VARIABLES

mkdir(info.save_path)
cd(info.save_path) 
save('alignment_variables','sound_info','sounds_per_file','sound_st','new_sound_st', 'sound_trials', 'sound_condition_array','info','passive_frames','conditions_per_speaker');
save('passive_frames','passive_frames')


%align and make imaging st

before_frames = 6;
after_frames = 91;
info.server = {'W:'}; %/Volumes/Runyan5
info.mouse_date = {'KN8-3L/2026-06-26'};

[imaging_st,temp] = align_passive_imagingst_updated_noise(info,before_frames,after_frames);
