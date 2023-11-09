% [sound_outputs,trialConditions, condition_onset_array_all]=
function find_spkr_output_task_new(server,mouse,date,alignment_info,spkr_channel_number,string,distance_between_sounds,threshold_spk,distance_within_sounds) 
%pc=1 if windows, any other number if mac
cd(strcat(server,'\Connie\RawData\',num2str(mouse),'\wavesurfer\',num2str(date)));
sync_dir = dir(strcat('*',string,'*.abf'));
num_files = length(sync_dir);
file_ind = 0;

for file = 1:num_files
    file_ind = file_ind +1
    [sync_data,sync_sampling_interval,~]  = abfload(sync_dir(file).name);
    sync_sampling_rate = 1/sync_sampling_interval*1e6;
    sync_data=sync_data';
    sync_data=double(sync_data);
    rawSounds=sync_data(spkr_channel_number,:);
    frames_times{file} = alignment_info(file).frame_times; 
    x = rescale(rawSounds(1,:),-1,1);
    [upperenv lowerenv] = envelope(x, 'linear');% problem
    value_n = 5;
    xmin = ordfilt2(x, 1, true(value_n));
    xmax = ordfilt2(x, value_n*value_n, true(value_n));
    figure(1);clf; 
    hold on;plot(x,LineWidth=0.05);plot(upperenv);plot(lowerenv); hold off
    value = 5;
    x1 = tanh(value*(upperenv-lowerenv)); %scale amplitude nonlinearly
    %x1 = tanh(value*(xmax-xmin)); %scale amplitude nonlinearly

    value2 = 40;
    x2 = conv2(x1,value2);%gaussian filter
    threshold = 1;%30;
    x3 = heaviside(threshold-x2);
    figure(2);clf ;hold on; plot(x2);plot(x3); hold off
    x4 = x3.*x;
    figure(3);clf; plot(x4)
    derivative = diff(x3);
    ups = find(derivative>0);
    down = find(derivative<0);

end
