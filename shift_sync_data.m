function shift_sync_data(data,file_trial_ids,digidata_its,file_estimated_trial_info)
[trial_its,trial_its_time] = virmen_it_rough_estimation(data); 

for file = 1:length(digidata_its)
%determine the timing of each virmen iteration in the session
iterations_in_time = data.data(1,:).*(86400);
iterations_in_time = iterations_in_time.*digidata_its(file).sync_sampling_rate;
it_ids = 1:length(data.data(1,:));

start_trial_number = file_trial_ids(file,1); 
end_trial_number = file_trial_ids(file,2);

%get the iterations that are within the imaging frames
first_file_it = digidata_its(file).locs(find(digidata_its(file).locs == file_estimated_trial_info(file).start_trials_digidata_time(file_trial_ids(file,3)))); %first starting it within imaging frames
last_file_it = digidata_its(file).locs(find(digidata_its(file).locs == file_estimated_trial_info(file).start_trials_digidata_time(file_trial_ids(file,4)))); %first starting it within imaging frames
%possible_iterations = digidata_its(file).locs(first_file_it): digidata_its(file).locs(last_file_it);


% case 1- I have a positive iteration! use this to determine shift and to
% number iterations accordingly
if ~isempty(digidata_its(file).pos_loc)
    pos_peak_id = digidata_its(file).pos_pks*10e4;
    pos_peak_pos = digidata_its(file).pos_loc(1); %use first one for now
    if pos_peak_id(1) == 1000000 %very first file has large positive peak
        shift = pos_peak_pos; 
        possible_it_times = iterations_in_time+shift-iterations_in_time(1); %assign iterations times with added shift
        possible_iterations = trial_its.start_trial_its(start_trial_number):trial_its.end_iti_its(end_trial_number); %limit iterations to ones within imaging frames
        %test distance of unfinished sounds and last iteration of trial to
        %see if they are close together otherwise add a shift of 1
        %iteration
    else
    
    
else
    % first say that the first iteration time in the first full trial is in the
    % position of the first iteration in the trua dataset (virmen data)
    possible_digidata_time_start = start_trials_digidata_time(2);%(start_trial_number); 
    %it_time_gaps
    %use time difference and assume first time difference to be exactly zero
    it_times_this_file = [data.data(1,possible_iterations)-data.data(1,possible_iterations(1))].*86400;% difference with first iteration for each iteration in seconds!
    it_times_this_file = it_times_this_file*sync_sampling_rate;
    
    %start by assuming that the estimated digidata time alings perfectly with
    %true data- here is where I figure out how to shift over time
    possible_it_times = it_times_this_file + possible_digidata_time_start;
end
%start_trial_its is the iti #
%possible iterations within file
possible_iterations = start_trial_its(start_trial_number): end_iti_its(end_trial_number);



ex_data = abfload(strcat(digidata_its(file).directory));
figure(55);clf;
title(strcat('Detected trial events file # ', num2str(file)));
hold on;aa = plot(ex_data(:,task_info.channel_number(1)));bb = plot(ex_data(:,task_info.channel_number(2)),'color',[0.7 0.7 0.7]);  cc = plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-b');dd = plot(rescale(ex_data(:,task_info.channel_number(4)),-1,0),'-m');
a = plot(end_trials_digidata_time,0,'*c');b = plot(start_trials_digidata_time,0,'*g');c= plot(end_iti_digidata_time,0,'*y'); movegui(gcf,'center');
%plot(rescale(ex_data(:,task_info.channel_number(3)),-1,0),'-r');
legend([aa bb cc dd  a(1) b(1) c(1)],'Imaging frames','Virmen its','Speaker 1','Speaker 2', 'end trial', 'start trial', 'end iti')
if length(task_info.channel_number)>4
    plot(rescale(ex_data(:,task_info.channel_number(5)),-1,0),'-r')
end
hold off;


end
%figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times,0,'*c');hold off; movegui(gcf,'center');
figure(); hold on;plot(ex_data(:,6)); plot(rescale(ex_data(:,4),-1,0));plot(possible_it_times(possible_iterations),0,'*c');hold off; movegui(gcf,'center');
end