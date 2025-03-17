%new bootstrap
function [sound_onsets_all,alignment_frames_all] = get_sound_trials_onsets(info,savepath)


for m = 1:length(info.mouse_date)
m
ss = info.serverid(m);
ss = ss {1,1};
base_path = strcat(num2str(ss),'\Connie\ProcessedData\',num2str(info.mouse_date{1,m}),'\');

%1) load data! 
load([base_path 'dff.mat']);
deconv = load([base_path 'deconv\deconv.mat']);
deconv = deconv.deconv;
load([base_path 'context_stim\60\context_tr.mat']);
load([base_path 'context_stim\60\bad_frames.mat']);
load([base_path 'passive\passive_frames.mat']); %load([base_path 'passive\passive_frames.mat']);


%check to make sure there are 3 repeats for each trial
rem_3 = rem(length(passive_frames.trial_num),passive_frames.trial_num(end));
if rem_3 > 0
    fprintf(['Dataset ', num2str(m),' has excluded trials\n']); %\n next line
%     keyboard
    % find sounds not in repeats of 3
    temp = []; for n = 1:length(unique(passive_frames.trial_num)); a = find(passive_frames.trial_num == n); temp = [temp, length(a)];end
    [aa,bb] = unique(passive_frames.trial_num);
    not_3 = bb(find(diff(temp)<0))+3; %sounds not in repeats of 3 - assumes previous sounds are a repeat of 3
    excluded_trials = passive_frames.trial_num(not_3);

    if length(not_3) < 1 %to deal with very first trial!
        excluded_trials = 1;
    end
else
    excluded_trials = [];
end

if info.data_type_boot == 1
    data = dff;
else
    data = deconv;
end

% get the trials that don't have stim
context_num = 2;
[control_output,opto_output] = find_control_frames(passive_frames,bad_frames,context_tr,context_num); %last is context number
alignment_based_on_opto = [bad_frames(context_tr{context_num,2},1),bad_frames(context_tr{context_num,2},2)]; %last number is ctrl (2) and stim (1)

%find trials without any stim!
[repeatloc,first_repeat] = unique(passive_frames.trial_num);

first_repeat = first_repeat(setdiff(repeatloc,excluded_trials));
sound_only_output = setdiff(first_repeat,[control_output;opto_output]); %frames that dont have opto or control just the sound!

%adjust alignment frames (based on bad_frames so I dont make adjustments
%using control and sound only trials
% sound_onset_frame = passive_frames.corr_frames([control_output;sound_only_output],1);
% alignment_frames = [sound_onset_frame,sound_onset_frame-2];
sound_onsets = [control_output;sound_only_output];
%alignment_frames = [passive_frames.corr_frames(:,1),passive_frames.corr_frames(:,1)-2]; %gets rid of any interpolation (if using -1 +2)
alignment_frames = [passive_frames.corr_frames(:,1)-1,passive_frames.corr_frames(:,1)+2]; %gives close interpolation as opto (if using -1 +2 for alignment)
alignment_based_on_control_opto = [bad_frames(context_tr{context_num,2},1),bad_frames(context_tr{context_num,2},2)];

%replace alignment sound frames with control opto frames for identical
%alignment!
alignment_frames(control_output,:) = alignment_based_on_control_opto;

%use bad_frames to make sure alignment is identical!


%perform for each location!

    sound_onsets_all{m} = sound_onsets;
    alignment_frames_all{m} = alignment_frames;
    control_output_all{m} = control_output;

end
%save variables
if ~isempty(savepath)
    mkdir([savepath '\mod'])
    cd([savepath '\mod'])

    save('sound_onsets_all','sound_onsets_all');
    save('alignment_frames_all','alignment_frames_all');
    save('control_output_all','control_output_all');
end