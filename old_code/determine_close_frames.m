function file_estimated_trial_info = determine_close_frames (alignment_info,file_estimated_trial_info)
%assumes structure field names are 1:4
% 'start_iti_digidata_timestart_trials_digidata_timeend_trials_digidata_timeend_iti_digidata_time'

field_names = fieldnames(file_estimated_trial_info(1));
for file = 1:length(file_estimated_trial_info)
    for te = 1:4
        trial_event = getfield(file_estimated_trial_info(file),field_names{te});
        for t = 1:length(trial_event)
            [val,closest_frame] = min(abs(trial_event (t) - alignment_info(file).frames_times));
            if te == 1
                if val < .02*alignment_info(file).sync_sampling_rate %FIND IF THERE IS FRAME 200 sec or .2ms from this one
                    file_estimated_trial_info(file).frame_start_iti(t) = closest_frame;
                else
                    file_estimated_trial_info(file).frame_start_iti(t) = nan;
                end
            elseif te == 2
                if val < .02*alignment_info(file).sync_sampling_rate %FIND IF THERE IS FRAME 200 sec or .2ms from this one
                file_estimated_trial_info(file).frame_start_trial(t) = closest_frame;
            else
                file_estimated_trial_info(file).frame_start_trial (t)= nan;
                end
            elseif te ==3
                if val < .02*alignment_info(file).sync_sampling_rate %FIND IF THERE IS FRAME 200 sec or .2ms from this one
                file_estimated_trial_info(file).frame_end_trial(t) = closest_frame;
            else
                file_estimated_trial_info(file).frame_end_trial(t) = nan;
                end
            elseif te ==4
                if val < .02*alignment_info(file).sync_sampling_rate %FIND IF THERE IS FRAME 200 sec or .2ms from this one
                file_estimated_trial_info(file).frame_end_iti(t) = closest_frame;
            else
                file_estimated_trial_info(file).frame_end_iti(t) = nan;
                end
            end

        end
    end
end
