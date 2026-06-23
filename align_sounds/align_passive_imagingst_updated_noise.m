function [imaging_st,temp] = align_passive_imagingst_updated_noise(info,before_frames,after_frames)
%updated to make sure sounds are aligned with bad_frames (sometimes
%slightly after)
% % before_frames = 7; %same as task
% % after_frames = 90;

imaging_st ={};
% inputs: info (mouse dates), number of frames in front, whether or not to
% get rid of interpolation during stim
for m = 1:length(info.mouse_date)
    m
    ss = info.server(m);
    ss = ss {1,1};
    base_path = strcat(num2str(ss),'\Connie\ProcessedData\',num2str(info.mouse_date{1,m}),'\');
    passive_savepath = strcat(num2str(ss),'\Connie\ProcessedData\',num2str(info.mouse_date{1,m}),'\passive\');

    
    %1) load data! 
    load([base_path 'dff.mat']);
    load([base_path 'z_dff.mat']);
    
    deconv = load([base_path 'deconv\deconv.mat']);
    deconv = deconv.deconv;
%     load([base_path 'context_stim\60\context_tr.mat']);
%     load([base_path 'context_stim\60\bad_frames.mat']);
bad_frames = [];

    load([base_path 'passive\passive_frames.mat']);
    load([base_path 'corrected_velocity.mat']); %contains 3 channels
%     load([base_path 'velocity_vector.mat']); %combines channels 1 & 2

%     oldFolderName = passive_savepath; % Full path to the current folder
%     newFolderName = strcat(passive_savepath, 'original'); % Generate a new name
%     
%     % Rename the folder
%     movefile(oldFolderName, newFolderName);


    %get rid of NAN trials
%     nan_trials = find(isnan(passive_frames.corr_frames(:,2)));
%     passive_frames.trial_num(nan_trials) = [];
%     nan_trials = find(isnan(passive_frames.corr_frames(:,1)));
%     passive_frames.trial_num(nan_trials) = [];
    
    %check to make sure there are 3 repeats for each trial
    excluded_trials = [];
    
    
    % get the trials that don't have stim
%     [control_output,opto_output] = find_control_frames(passive_frames,bad_frames,context_tr,2); %last is context number
    
    %find trials without any stim!
    [repeatloc,first_repeat] = unique(passive_frames.trial_num);
    keep_trials = setdiff(repeatloc, excluded_trials);
    keep_idx = ismember(repeatloc, keep_trials);
    
    first_repeat = first_repeat(keep_idx);
    repeatloc = repeatloc(keep_idx);

    [repeatloc, first_repeat] = unique(passive_frames.trial_num);

    %in case two trial types are identical and next to each other
    extra_repeatloc = [];
    extra_first_repeat = [];
    
    for i = 1:length(repeatloc)
    
        idx = find(passive_frames.trial_num == repeatloc(i));
    
        % first occurrence already stored by unique
        if length(idx) > 3
    
            % add every occurrence after the third
            extra_first_repeat = [extra_first_repeat idx(4:end)];
            extra_repeatloc = [extra_repeatloc repmat(repeatloc(i),1,length(idx)-3)];
    
        end
    end
    
    repeatloc = [repeatloc; extra_repeatloc(:)];
    first_repeat = [first_repeat; extra_first_repeat(:)];
    
    [first_repeat,sort_idx] = sort(first_repeat);
    repeatloc = repeatloc(sort_idx);
        
%     first_repeat = first_repeat(setdiff(repeatloc,excluded_trials));
%     sound_only_output = setdiff(first_repeat,[control_output;opto_output]); %frames that dont have opto or control just the sound!
    
    %adjust alignment frames (based on bad_frames so I dont make adjustments
    %using control and sound only trials
%     sound_onsets = [control_output;sound_only_output];
%     alignment_frames = [passive_frames.corr_frames(:,1),passive_frames.corr_frames(:,1)-2]; %gets rid of any gap (if using -1 +2)
%     alignment_frames = [passive_frames.corr_frames(:,1),passive_frames.corr_frames(:,1)+2]; %gives same gap as opto (if using -1 +2) 
    alignment_frames = [passive_frames.corr_frames(:,1),passive_frames.corr_frames(:,1)]; %gives same gap as opto (if using -1 +2) 
    
    
    % % x1=frames_before_event+1;
    % % x2=1;
    % % y1=2;
    % % y2=frames_after_event+2;
    % % bfint= x-x1:x-x2; 
    % %    afint= y+y1:y+y2;
    
    context_num = 2; %1 task, 2 passive, 3 spont
    
    temp = [];
    
    %to keep track of files take difference of original frame and find when
    %the difference in negative (add 1)
    file_starts = [find(diff(passive_frames.frames(:,1))<0)+1];
    file_starts = [file_starts;length(passive_frames.frames)];

    %update passive frames to make sure they are aligned with photostim
    passive_frames.corr_frames_updated = passive_frames.corr_frames;

    %decide how to adjust based on photostim
    to_add = ones(length(first_repeat),1);
    to_add = to_add *0;
    for trial = 1:length(first_repeat)
        if ~isempty(bad_frames) %determine stim trials and control?
                if any(find(ismember(passive_frames.corr_frames(first_repeat(trial),1)-1,bad_frames(context_tr{context_num,1}))) ==1) 
                    passive_frames.corr_frames_updated(first_repeat(trial),1) = passive_frames.corr_frames(first_repeat(trial),1)-1;
                    to_add(trial) = [1];
                elseif any(find(ismember(passive_frames.corr_frames(first_repeat(trial),1)-2,bad_frames(context_tr{context_num,1}))) ==1)
                    passive_frames.corr_frames_updated(first_repeat(trial),1) = passive_frames.corr_frames(first_repeat(trial),1)-2;
                    to_add(trial) = [2];
                end
    
                if any(find(ismember(passive_frames.corr_frames(first_repeat(trial),1)-1,bad_frames(context_tr{context_num,2}))) ==1) 
                    passive_frames.corr_frames_updated(first_repeat(trial),1) = passive_frames.corr_frames(first_repeat(trial),1)-1;
                    to_add(trial) = [1];
                elseif any(find(ismember(passive_frames.corr_frames(first_repeat(trial),1)-2,bad_frames(context_tr{context_num,2}))) ==1)
                    passive_frames.corr_frames_updated(first_repeat(trial),1) = passive_frames.corr_frames(first_repeat(trial),1)-2;
                    to_add(trial) = [2];
                end
        end 
    end
    
    %initialize structure
    imaging = {};trial_lengths = []; temp = [];
    for trial = 1:length(first_repeat)
        trial
        trial_id = repeatloc(trial);

        % Find all rows corresponding to this trial
        trial_rows = find(passive_frames.trial_num == trial_id);
    
        start_trial_frame = ...
            passive_frames.corr_frames_updated(trial_rows(1),1) - before_frames;
    
        end_trial_frame = ...
            passive_frames.corr_frames_updated(trial_rows(end),2) + after_frames;
        
        trial_length = length(start_trial_frame:end_trial_frame);
        trial_lengths = [trial_lengths,trial_length];
    % for each trial in this folder/acquisition find the start and end iteration #
        imaging(trial).start_it = [];
        imaging(trial).end_it =  [];
        imaging(trial).iti_start_it = [];
        imaging(trial).iti_end_it =  [];
    
            %virmen trial info! for all trials even without imaging
        virmen_trial_info(trial).correct = 0; %correct or not for that trial
        virmen_trial_info(trial).left_turn = 0; %left turn or not for that trial
    
        virmen_trial_info(trial).condition = passive_frames.condition(first_repeat(trial)); %stimulus condition for that trial
        if ~isempty(bad_frames) %determine stim trials and control?
            if any(find(ismember(passive_frames.corr_frames_updated(first_repeat(trial),1)-1,bad_frames(context_tr{context_num,1}))) ==1) ||  any(find(ismember(passive_frames.corr_frames_updated(first_repeat(trial),1),bad_frames(context_tr{context_num,1}))) ==1)|| any(find(ismember(passive_frames.corr_frames_updated(first_repeat(trial),1)-2,bad_frames(context_tr{context_num,1}))) ==1)
                virmen_trial_info(trial).is_stim_trial = 1;
                temp = [temp,trial];
            else
                virmen_trial_info(trial).is_stim_trial = 0;
            end
        end 
    
        imaging(trial).virmen_trial_info = virmen_trial_info(trial);
    
        %zero for everything except running velocity and stimulus!
        movement_in_virmen_time(trial).y_position = [];%zeros(1,trial_length);
        movement_in_virmen_time(trial).x_position = [];%zeros(1,trial_length);
        movement_in_virmen_time(trial).view_angle = [];%zeros(1,trial_length);
        movement_in_virmen_time(trial).x_velocity = [];%corrected_velocity(2,start_trial_frame:end_trial_frame); %2 is roll
        movement_in_virmen_time(trial).y_velocity = [];%corrected_velocity(1,start_trial_frame:end_trial_frame); %1 is pitch
        movement_in_virmen_time(trial).cw = [];%passive_frames.condition(first_repeat(trial)); %current world
        movement_in_virmen_time(trial).is_reward = [];%zeros(1,trial_length);
        movement_in_virmen_time(trial).in_ITI = [];%zeros(1,trial_length); %add a one or 2 second ITI after sounds are done played
        
        if ~isnan(start_trial_frame) && ~isnan(end_trial_frame)
            movement_in_imaging_time.y_position = zeros(1,trial_length);
            movement_in_imaging_time.x_position = zeros(1,trial_length);
            movement_in_imaging_time.y_velocity = corrected_velocity(2,start_trial_frame:end_trial_frame); %2 is roll
            movement_in_imaging_time.x_velocity = corrected_velocity(1,start_trial_frame:end_trial_frame); %1 is pitch
            movement_in_imaging_time.view_angle = zeros(1,trial_length);
            movement_in_imaging_time.is_reward = zeros(1,trial_length);
            movement_in_imaging_time.in_ITI = zeros(1,trial_length);
            movement_in_imaging_time.pure_tones =zeros(1,trial_length);
            movement_in_imaging_time.maze_frames = 1:(passive_frames.corr_frames_updated(first_repeat(trial)+2,2)+1-start_trial_frame);
            movement_in_imaging_time.iti_frames = (passive_frames.corr_frames_updated(first_repeat(trial)+2,2)+1-start_trial_frame)+1:trial_length;
            movement_in_imaging_time.reward_frames = zeros(1,trial_length);
    
            %binarize stimulus onsets
            temp_stimulus = zeros(1,trial_length);
    %         stimulus_onsets = [passive_frames.corr_frames(first_repeat(trial),1)-start_trial_frame+1:passive_frames.corr_frames(first_repeat(trial),2)-start_trial_frame+1,passive_frames.corr_frames(first_repeat(trial)+1,1)-start_trial_frame+1:passive_frames.corr_frames(first_repeat(trial)+1,2)-start_trial_frame+1,passive_frames.corr_frames(first_repeat(trial)+2,1)-start_trial_frame+1:passive_frames.corr_frames(first_repeat(trial)+2,2)-start_trial_frame+1];
    
%             stimulus_onsets = [passive_frames.corr_frames_updated(first_repeat(trial),1)-start_trial_frame+1:passive_frames.corr_frames_updated(first_repeat(trial),2)-start_trial_frame+1,passive_frames.corr_frames_updated(first_repeat(trial)+1,1)-start_trial_frame+1:passive_frames.corr_frames_updated(first_repeat(trial)+1,2)-start_trial_frame+1,passive_frames.corr_frames_updated(first_repeat(trial)+2,1)-start_trial_frame+1:passive_frames.corr_frames_updated(first_repeat(trial)+2,2)-start_trial_frame+1];
            stimulus_onsets = [];
            
            for r = 1:length(trial_rows)
            
                stimulus_onsets = [stimulus_onsets, ...
                    passive_frames.corr_frames_updated(trial_rows(r),1)-start_trial_frame+1 : ...
                    passive_frames.corr_frames_updated(trial_rows(r),2)-start_trial_frame+1];
            
            end      

            stimulus_onsets = stimulus_onsets(~isnan(stimulus_onsets));
            temp_stimulus(stimulus_onsets) = 1;
            
            movement_in_imaging_time.stimulus = temp_stimulus;
        
            imaging(trial).movement_in_virmen_time = movement_in_virmen_time(trial);
        
            imaging(trial).frame_id = start_trial_frame:end_trial_frame;% virmen iterations that fit within this trial including ITI
        
            maze_start_frame = start_trial_frame;
            iti_end_frame = end_trial_frame;
        
            %finding the closest frames per trial epoch
            imaging(trial).frame_id_events.maze =start_trial_frame:passive_frames.corr_frames_updated(first_repeat(trial)+2,2);
            imaging(trial).frame_id_events.iti = passive_frames.corr_frames_updated(first_repeat(trial)+2,2)+1:passive_frames.corr_frames_updated(first_repeat(trial)+2,2)+after_frames;
        
            imaging(trial).dff = dff(:,start_trial_frame:end_trial_frame);
            imaging(trial).z_dff = z_dff(:,start_trial_frame:end_trial_frame);
            imaging(trial).deconv = deconv(:,start_trial_frame:end_trial_frame);
            
        %     imaging(trial).relative_frames = maze_start_frame+previous_frames_sum:iti_end_frame+previous_frames_sum;
            
            imaging(trial).file_num = find(first_repeat(trial) < file_starts,1,'first');
        
            movement_in_imaging_time.turn_frame = 0; 
            %save movement in imaging frames
            imaging(trial).movement_in_imaging_time = movement_in_imaging_time; 
            
            %good trial: the entire trial is usable
            imaging(trial).good_trial = 1;
        else
            imaging(trial).good_trial = 0;
        end


    end
    to_add
%save variables
if ~isempty(passive_savepath)
    cd(passive_savepath)
    save('imaging','imaging');  
    save('frames_used',"before_frames","after_frames");
end

imaging_st{1,m} = imaging;
end
