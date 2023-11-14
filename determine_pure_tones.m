function [condition_onset_array_all,pure_tones_trial,correct_pure_tones,incorrect_pure_tones] = determine_pure_tones(pure_tone_signal,sync_sampling_rate,correct_length,incorrect_length,condition_onset_array_all)
pos_peaks = []; pos_peaks = [];
for s = 1:size(pure_tone_signal,1)
    [~,pos_peaks2] = findpeaks(abs(diff(pure_tone_signal(s,:))),'MinPeakHeight', 0.005);
    pos_peaks = [pos_peaks,pos_peaks2];
end
correct_length = correct_length*sync_sampling_rate;
incorrect_length = incorrect_length*sync_sampling_rate;
time_after = [0.1*sync_sampling_rate,1*sync_sampling_rate]; %examples were 0.188 0.192 0.193 (incorrect), .523 .664 0.687 correct
target_differences = [correct_length:correct_length*.02+correct_length];
target_differences_incorrect = [incorrect_length-incorrect_length*.02:incorrect_length];

combos = nchoosek(1:length(pos_peaks),2);
correct_pure_tones = [];incorrect_pure_tones = [];pure_tones_trial = [];trial_options=[];trial_options2=[]; possible_values = [];
all_differences = abs(pos_peaks(combos(:,2)) - pos_peaks(combos(:,1)));
in_diff = find(all_differences>target_differences_incorrect(1) & all_differences < target_differences_incorrect(end));
cor_diff = find(all_differences>target_differences(1) & all_differences < target_differences(end));

good_diff= [ in_diff,cor_diff];
for t = 1:size(condition_onset_array_all.VR_sounds,1)
    t
    for ii = 1:length(good_diff)
    
    possible_values = [];

            %find differences between each peak value
            i = good_diff(ii);
            difference = all_differences(i);
            % check if difference is equal to target differences
            if any(difference == target_differences_incorrect)&& min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)<time_after(2) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)>time_after(1)
                %incorrect_pure_tones = [incorrect_pure_tones;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
%                 if size(pure_tones_trial,1)==(t) %deleting correct tone if I find incorrect (since correct would be inside incorrect)
%                     correct_pure_tones(t,:) = [0,0];
%                 end
                trial_options = [trial_options;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
                pure_tones_trial(t,:) = [0,trial_options(1,:),t];
                 
                incorrect_pure_tones(t,:) = [trial_options(1,:)];
            
            elseif any(difference == target_differences) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)<time_after(2) && min([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))]) - condition_onset_array_all.VR_sounds(t,3)>time_after(1) && size(pure_tones_trial,1)==(t-1)
                %correct_pure_tones = [correct_pure_tones;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
                trial_options2 = [trial_options2;sort([pos_peaks(combos(i,1)), pos_peaks(combos(i,2))])];
                pure_tones_trial(t,:) = [1,trial_options2(1,:),t];
                correct_pure_tones(t,:) = [trial_options2(1,:)];
            %end
            end
            
    end
    trial_options=[];
    trial_options2=[];
    if size(pure_tones_trial,1) < t
       pure_tones_trial(t,:) = [nan,nan,nan,t];
    end
    condition_onset_array_all.ITI_sounds(t,:) = pure_tones_trial(t,:);
end

figure(112);clf; hold on; title('ITI Pure tones with NaN')
hold on
plot(diff(pure_tone_signal(1,:)));
plot((pure_tone_signal(1,:)));
if length(find(isnan(pure_tones_trial(:,1))))>0  
    plot(condition_onset_array_all.VR_sounds([find(isnan(pure_tones_trial(:,1)))],3),0,'*m');
end
hold off
legend('Pure Tones Diff','Pure Tones regular','Last sound near NaN')

% plot(correct_pure_tones(:,1),0,'*b');
% plot(correct_pure_tones(:,2),0,'+b'); 
% plot(incorrect_pure_tones(:,1),0,'*r');
% plot(incorrect_pure_tones(:,2),0,'+r');
% legend('','Correct','','Incorrect','')
% hold off