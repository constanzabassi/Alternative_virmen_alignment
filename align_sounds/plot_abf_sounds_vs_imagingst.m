%check imaging st vs actual sounds/opto
%load necessary files
mouse_ex = ['HA1-00'];%'HA11-1R';
date_ex = '2023-06-29';%%'2023-05-05';
server = 'V';
og = load([server ':\Connie\ProcessedData\' mouse_ex '\' date_ex '\passive\original\imaging.mat'])
load([server ':\Connie\ProcessedData\' mouse_ex '\' date_ex '\alignment_info.mat'])
load([server ':\Connie\ProcessedData\' mouse_ex '\' date_ex '\bad_frames.mat'])

%%
frame_lengths = cellfun(@(x) length(x),{alignment_info.frame_times});
total_lengths = [0, cumsum(frame_lengths)];
file_num = 11;%9; %;11;
m=5;%20;

to_sub = total_lengths(file_num);
[sync,a] = abfload([server ':\Connie\RawData\' mouse_ex '\wavesurfer\' date_ex '\' alignment_info(file_num).sync_id]);
opto_chan = 5;
trials = [4,7,10,20];

figure(1);clf;
for i = 1:4

frame_range = og.imaging(trials(i)).frame_id(1)-to_sub:og.imaging(trials(i)).frame_id(end)-to_sub(end);
sync_times = alignment_info(file_num).frame_times(frame_range);

subplot(2,2,i);
hold on
sound_id = og.imaging(trials(i)).virmen_trial_info.condition;
if sound_id ==2
    chan = 8;
else 
    chan = 4;
end
plot(rescale(sync(sync_times,chan),0,1),'-b');
plot(rescale(sync(sync_times,opto_chan),0,1),'-r');
plot(og.imaging(trials(i)).movement_in_imaging_time.stimulus,'--k');
%xlim([1 100])
hold off
end

figure(2);clf;
for i = 1:4

frame_range = imaging_st{1,m}(trials(i)).frame_id(1)-to_sub:imaging_st{1,m}(trials(i)).frame_id(end)-to_sub(end);
sync_times = alignment_info(file_num).frame_times(frame_range);

subplot(2,2,i);
hold on
sound_id = imaging_st{1,m}(trials(i)).virmen_trial_info.condition;
if sound_id ==2
    chan = 8;
else 
    chan = 4;
end
plot(rescale(sync(sync_times,chan),0,1),'-b');
plot(rescale(sync(sync_times,opto_chan),0,1),'-r');
plot(imaging_st{1,m}(trials(i)).movement_in_imaging_time.stimulus,'--k');
%xlim([1 100])
hold off
end

figure(3);clf;

for i = 1:4

frame_range1 = og.imaging(trials(i)).frame_id(1)-to_sub:og.imaging(trials(i)).frame_id(end)-to_sub(end);
sync_times1 = alignment_info(file_num).frame_times(frame_range1);

frame_range = imaging_st{1,m}(trials(i)).frame_id(1)-to_sub:imaging_st{1,m}(trials(i)).frame_id(end)-to_sub(end);
sync_times = alignment_info(file_num).frame_times(frame_range);

subplot(2,2,i);
hold on
sound_id = imaging_st{1,m}(trials(i)).virmen_trial_info.condition;
if sound_id ==2
    chan = 8;
else 
    chan = 4;
end
plot(rescale(sync(sync_times,chan),0,1),'-b');
plot(rescale(sync(sync_times,opto_chan),0,1),'-r');
plot(rescale(sync(sync_times1,chan),0,1),'-c');
plot(rescale(sync(sync_times1,opto_chan),0,1),'-y');

plot(og.imaging(trials(i)).movement_in_imaging_time.stimulus,'-g');
plot(imaging_st{1,m}(trials(i)).movement_in_imaging_time.stimulus,'-m');
%xlim([1 100])
hold off
end