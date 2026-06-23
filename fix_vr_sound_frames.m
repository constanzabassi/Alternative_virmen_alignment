function vr_sound_frames_updated = fix_vr_sound_frames(vr_sound_frames, imaging, alignment_info)

frame_lengths = [];
frame_lengths = cellfun(@length,{alignment_info.frame_times}); %across all imaged files 
frame_lengths = [0,cumsum(frame_lengths)];

vr_sound_frames_updated = vr_sound_frames;

n_sounds = size(vr_sound_frames.frames,1);

% ------------------------------------------------------------
% 1. Identify valid imaging trials and build mapping
% ------------------------------------------------------------
valid_imaging_idx = find(~cellfun(@isempty, {imaging.good_trial}));

% mapping: imaging index -> consecutive trial id
imaging_to_newtrial = nan(1, length(imaging));
imaging_to_newtrial(valid_imaging_idx) = 1:length(valid_imaging_idx);

% ------------------------------------------------------------
% 2. Assign each sound to an imaging trial via frame overlap
% ------------------------------------------------------------
new_trial_num = nan(n_sounds,1);

for i = 1:length(valid_imaging_idx)

    tr = valid_imaging_idx(i);

    if isempty(imaging(tr).frame_id)
        continue
    end

    tr_start = imaging(tr).frame_id(1) + frame_lengths(imaging(tr).file_num);
    tr_end   = imaging(tr).frame_id(end)+ frame_lengths(imaging(tr).file_num);

    sound_idx = vr_sound_frames.corr_frames(:,1) >= tr_start & ...
                vr_sound_frames.corr_frames(:,2) <= tr_end;

    if any(sound_idx)

        % sanity check: all sounds should come from same original VR trial
        old_trials = unique(vr_sound_frames.trial_num(sound_idx));
        if numel(old_trials) > 1
            warning('Imaging trial %d spans multiple VR trial_nums', tr);
        end

        new_trial_num(sound_idx) = valid_imaging_idx(imaging_to_newtrial(tr));
    end
end

% ------------------------------------------------------------
% 3. Remove sounds not assigned to imaging trials
% ------------------------------------------------------------
keep = ~isnan(new_trial_num);

fields = fieldnames(vr_sound_frames);

for f = 1:length(fields)
    data = vr_sound_frames.(fields{f});
    if size(data,1) == n_sounds
        vr_sound_frames_updated.(fields{f}) = data(keep,:);
    end
end

% ------------------------------------------------------------
% 4. Update trial numbering
% ------------------------------------------------------------
vr_sound_frames_updated.trial_num = new_trial_num(keep); %relative to virmen

end
