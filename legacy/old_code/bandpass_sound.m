function sound_peaks = bandpass_sound(sound_matrix,freq_range,min_distance_between,min_peak_height,sync_sampling_rate)
for ss = 1:size(sound_matrix,1)
    filtered_sound(ss,:) = bandpass(sound_matrix(ss,:),freq_range,sync_sampling_rate); %[2 4] [0.2 2]
end

sound_peaks = [];

figure();
title(strcat('Bandpass signal with frequencies between ', num2str(freq_range), ' Hz'))
for ss = 1:size(sound_matrix,1)
    subplot(size(sound_matrix,1),1,ss)
    [~, z] = findpeaks(abs(filtered_sound(ss,:)),'MinPeakDistance',min_distance_between,'MinPeakHeight',min_peak_height); %,0.0007/0.00095
     hold on
     minn = min(filtered_sound(ss,:));
     maxx = max(filtered_sound(ss,:));
     plot(rescale(sound_matrix(ss,:),minn,maxx));  plot(filtered_sound(ss,:));plot(z,abs(filtered_sound(ss,z)),'*r');
     hold off
     sound_peaks = [sound_peaks,z];
end