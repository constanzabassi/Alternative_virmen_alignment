function [trial_its,trial_its_time] = virmen_it_rough_estimation(data)
% use virmen iterations to get rough estimate of which iterations are
% specific trial events based on stereotyped gaps between them!
% also find iteration where mouse crosses 50 in the y position!- assumes a
% threshold of 350- takes more than these # iterations to do a trial
threshold = 350; %lowest value for this mouse is ~454 in 500 maze- distance between sounds on different trials


%look at iti to determine isITI and also start of trial
start_iti_its = find(diff(data.data(9,:))>0);
start_iti_its = start_iti_its +1; 
end_iti_its = find(diff(data.data(9,:))<0);
start_trial_its = end_iti_its + 1;
start_trial_its = sort([start_trial_its,1]); %making sure I have an it for the very first trial
end_trial_its = start_iti_its - 1;


trial_its.end_trial_its = end_trial_its;
trial_its.start_iti_its = start_iti_its;
trial_its.end_iti_its = end_iti_its;

%get complete trials
min_trials = min([length(end_trial_its),length(start_trial_its)]);
if length(start_trial_its) > min_trials
    start_trial_its(end) = [];
end
trial_its.start_trial_its = start_trial_its;

%% looking at iteration timing info! try to tie it back to the big sterotyped gaps that I see!
it_times = data.data(1,:).* 86400; %convert to seconds
it_time_gaps = diff(it_times);

%look at iti to determine isITI and also start of trial in iteration times
%did not use this really?
end_trial_its_time = it_times(find(it_time_gaps>.8));
start_iti_its_time = it_times(find(it_time_gaps>.8)+1);
incorrect_its_time = it_times(find(it_time_gaps >.8 & it_time_gaps < .95));
correct_its_its_time = it_times(find(it_time_gaps > .95));
temp_its = find(it_time_gaps > .25 & it_time_gaps < .51)+1;
start_trial_its_time = it_times([1,temp_its]); %adding one for very first trial

trial_its_time.start_trial_its_time = start_trial_its_time;
trial_its_time.end_trial_its_time = end_trial_its_time;
trial_its_time.start_iti_its_time = start_iti_its_time;
trial_its_time.correct_its_its_time = correct_its_its_time;
trial_its_time.incorrect_its_time = incorrect_its_time;

%look at sound onset first
temp_onset = sort([find(round(floor(data.data(3,:))) == 50),find(round(floor(data.data(3,:))) == 51),find(round(floor(data.data(3,:))) == 52)]);
%find ones that are too close together and delete them
temp = find(diff(temp_onset)< threshold);
temp_onset(temp+1) = [];
sound_trigger_its = temp_onset;
%make sure there is one trigger per trial! finding repeats getting rid of
%anything after first one
temp_sound_trigger = [];
for t = 1:min([length(trial_its.start_trial_its),length(trial_its.end_trial_its)])
    find(abs(trial_its.start_trial_its(t) - sound_trigger_its) == min(abs(trial_its.start_trial_its(t) - sound_trigger_its))) 
    temp_sound_trigger = [temp_sound_trigger, t];
end
[U, I] = unique(temp_sound_trigger, 'first');
x = 1:length(temp_sound_trigger);
x(I) = [];
sound_trigger_its(x) = [];% getting rid of repeats after first

trial_its.sound_trigger_its = sound_trigger_its;