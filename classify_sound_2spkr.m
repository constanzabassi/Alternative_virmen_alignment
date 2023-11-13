function [sound_struc, condition_array, onset_array, offset_array,classified_sounds] = classify_sound_2spkr (sound_signal,all_trial_sounds)

condition_array = []; onset_array = []; offset_array= []; tt = 0;sound_struc = {};
num_sounds = size(sound_signal,1);
classified_sounds = zeros(1,size(sound_signal,2));
  for t = 1:length(all_trial_sounds)
        current_sound = all_trial_sounds(t,1):all_trial_sounds(t,2);
        
        [~,spkr_loc] = max(sound_signal(:,all_trial_sounds(t,1)+1));
        
%         if any(norm_sounds(1,current_sound)>threshold_spk) && max(norm_sounds(2,current_sound))<threshold_spk %speaker 1
%             sounds(current_sound)=1;
%     
%         elseif any(norm_sounds(2,current_sound)>threshold_spk) && max(norm_sounds(1,current_sound))<threshold_spk %spekear 2
%             sounds(current_sound)=2;
%         end 
         if sum(sound_signal(setdiff(1:num_sounds,spkr_loc),all_trial_sounds(t,1)+2))<1 %making sure there are no signals in other speakers
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