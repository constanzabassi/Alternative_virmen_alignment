function [sound_struc, condition_array, onset_array, offset_array,classified_sounds] = classify_sound_2spkr (sound_signal,all_trial_sounds,mult_spkr)

condition_array = []; onset_array = []; offset_array= []; tt = 0;sound_struc = {};
num_sounds = size(sound_signal,1);
classified_sounds = zeros(1,size(sound_signal,2));
for t = 1:length(all_trial_sounds)
        current_sound = all_trial_sounds(t,1):all_trial_sounds(t,2);
        
        
        
%         if any(norm_sounds(1,current_sound)>threshold_spk) && max(norm_sounds(2,current_sound))<threshold_spk %speaker 1
%             sounds(current_sound)=1;
%     
%         elseif any(norm_sounds(2,current_sound)>threshold_spk) && max(norm_sounds(1,current_sound))<threshold_spk %spekear 2
%             sounds(current_sound)=2;
%         end 
if mult_spkr ==1
    spkr_loc = find(sound_signal(:,all_trial_sounds(t,1)+5));
    if length(spkr_loc) < size(sound_signal,1) && length(spkr_loc)>0 %sounds that are across ALL speakers are ITI!!
                     tt = tt+1;
                    sound_struc(tt).onset = all_trial_sounds(t,1);
                    sound_struc(tt).offset = all_trial_sounds(t,2);
                    sound_struc(tt).spkr_loc = spkr_loc;
                    if length(spkr_loc) > 1
                        sound_struc(tt).condition = 4;%sounds(all_trial_sounds(t,1)+2);
                    else
                        sound_struc(tt).condition = spkr_loc;
                    end
                    condition_array = [condition_array,sound_struc(tt).condition];
                    onset_array = [onset_array,all_trial_sounds(t,1)];
                    offset_array = [offset_array,all_trial_sounds(t,2)];
                    classified_sounds(current_sound) = sound_struc(tt).condition;
    end

else
    [~,spkr_loc] = max(sound_signal(:,all_trial_sounds(t,1)+1));
             if sum(sound_signal(setdiff(1:num_sounds,spkr_loc),all_trial_sounds(t,1)+2))<1 &&length(spkr_loc)>0 %making sure there are no signals in other speakers 
                 tt = tt+1;
                sound_struc(tt).onset = all_trial_sounds(t,1);
                sound_struc(tt).offset = all_trial_sounds(t,2);
                sound_struc(tt).condition = spkr_loc;%sounds(all_trial_sounds(t,1)+2);
                condition_array = [condition_array,sound_struc(tt).condition];
                onset_array = [onset_array,all_trial_sounds(t,1)];
                offset_array = [offset_array,all_trial_sounds(t,2)];
                classified_sounds(current_sound) = spkr_loc;
            end
end

end
