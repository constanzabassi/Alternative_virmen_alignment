function [all_sounds_together,new_sound_st,sounds_per_file] = binarize_passive_sounds(sound_st,sound_info,alignment_info,varargin)

%get dir ids for the alignment info!!
% Parse optional inputs
p = inputParser;
addParameter(p,'dir_string','passive',@(x) ischar(x) || isstring(x));
parse(p,varargin{:});

dir_string = char(p.Results.dir_string);

passive_temp_dir = cellfun(@(x) contains(x,dir_string),{alignment_info.sync_id},'UniformOutput',false);
passive_dir = find([passive_temp_dir{1,:}]); %gives the dir IDs of directory with passive in their name

%add up all frames in the prior alignment structures (so things line up
%with suite2p where all frames are concatenated together)
frame_lengths = cellfun(@length,{alignment_info.frame_times});

for file = 1:length(sound_st)

    if file > 1 || passive_dir(1) > 1
        previous_sum = sum(frame_lengths(1:passive_dir(file)-1)); %keep track of all previous frames
    else
        previous_sum = 0;
    end


    sounds_per_file(file).onsets = [sound_st(file).file.onset];%[[sound_st(file).file.onset],[sound_onsets{:,:}]];
    sounds_per_file(file).offsets = [sound_st(file).file.offset];%[[sound_st(file).file.offset],[sound_offsets{:,:}]];
    new_sound_st(file).condition = [sound_st(file).file.true_condition]';
    new_sound_st(file).og_trial_num = [sound_st(file).file.trial_num]';
    if file == 1
        new_sound_st(file).trial_num = [sound_st(file).file.trial_num]';
    else
        new_sound_st(file).trial_num = [sound_st(file).file.trial_num]'+sum(cellfun(@(x) sum(x(end)),{new_sound_st(file-1:-1:1).og_trial_num}));
    end
    
    
    binary_sounds = zeros(1,(max(sounds_per_file(file).offsets)+sound_info.sync_sampling_rate));
    binary_sounds_frames = [];
    min_distance = mode(diff(alignment_info(passive_dir(file)).frame_times))+(2*alignment_info(passive_dir(file)).sync_sampling_rate/1000); %minimum distance between frames in digidata time
    binary_sounds_imaging_time = zeros(1,length(alignment_info(passive_dir(file)).frame_times));
    
    for p = 1:min([length(sounds_per_file(file).offsets),length(sounds_per_file(file).onsets)])
        if p == 86
            a = 1;
        end
        %binarize signal based on digidata time
        binary_sounds(sounds_per_file(file).onsets(p):sounds_per_file(file).offsets(p)) = 1; %in digidata time
        %find onset and offset frames and binarize based on their positions
        binary_sounds_frames(p,:) =  [find_closest_frames(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).onsets(p),min_distance),find_closest_frames(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).offsets(p),min_distance)];
        %assign frames to sound_st
        new_sound_st(file).frames(p,:) =  [find_closest_frames(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).onsets(p),min_distance),find_closest_frames(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).offsets(p),min_distance)];
            
        if sound_info.distance_within_sounds < 1*sound_info.sync_sampling_rate
            %FIX FRAMES if they are at the very very start of the next sound
            %onset - if less than 2/3 of the sound then go to the next frame
            %%helps with sounds separated by 50ms!
            [~,dist] = find_closest_frames_adjusted(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).onsets(p),min_distance);
            [~,dist2] = find_closest_frames_adjusted(alignment_info(passive_dir(file)).frame_times,sounds_per_file(file).offsets(p),min_distance);
            if dist <round(mode(diff(alignment_info(passive_dir(file)).frame_times))*(2/3))
                new_sound_st(file).frames(p,1) = new_sound_st(file).frames(p,1)+1;
                binary_sounds_frames(p,1) = binary_sounds_frames(p,1)+1;
            end
            if dist2>round(mode(diff(alignment_info(passive_dir(file)).frame_times))*(2/3))
                new_sound_st(file).frames(p,2) = new_sound_st(file).frames(p,2)-1;
                binary_sounds_frames(p,2) = binary_sounds_frames(p,2)-1;
            end
        end
        new_sound_st(file).corr_frames(p,:) = new_sound_st(file).frames(p,:)+previous_sum;

        if all(~isnan(binary_sounds_frames(p,:)))
            binary_sounds_imaging_time(binary_sounds_frames(p,1):binary_sounds_frames(p,2)) = 1;
        end
    end
    figure(557);clf;
    hold on
    title(strcat('Frame times sounds in file # ', num2str(file)));
    plot(binary_sounds_imaging_time,'-k')
    hold off
    pause(2)
    
    sounds_per_file(file).binary_digidata_times = binary_sounds;
    sounds_per_file(file).binary_frame_times = binary_sounds_imaging_time;

end

all_sounds_together = [];
f = fieldnames(new_sound_st);
for i = 1:length(f)
    all_sounds_together.(f{i}) = cat(1,new_sound_st(1:length(sound_st)).(f{i}));%[new_sound_st(1:length(sound_st)).(f{i})];%{[new_sound_st(1:2).(f{i})]};
end

